#!/usr/bin/env python

import argparse
import xarray as xr
import logging
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy
import cartopy.crs as crs
import cartopy.feature as cfeature
from PIL import Image
from itertools import product

cartopy.config['data_dir'] = '/root/.local/share/cartopy'

plt.rcParams.update({'font.size': 8})

# Configure basic logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("PLOTTER")

# Time format
FMT_TIME = "%Y-%m-%d %H:%MZ"

class MapImage:
    def __init__(self,ncfile):
        self.key      = None
        self.dims     = None
        self.data     = None
        self.indexes  = None
        self.hires    = False
        ###
        ### Open netcdf
        ###
        self.ds = xr.open_dataset(ncfile)
        self.nt = self.ds.sizes['time']
        ###
        ### Plot config
        ###
        if self.ds.lon.cell_measures < 0.1 or self.ds.lat.cell_measures < 0.1:
            self.hires = True

    def set_key(self,key):
        if key in self.ds:
            self.key  = key
            self.dims = self.ds[key].dims
        else:
            logger.exception(f"Variable {key} not found")
            raise RuntimeError("Aborting program due to error: key not found")

    def key_indexes(self, **required_indexes):
        indexable_dims = ['time','lev','bin','fl','layer']
        indexes = {dim: [0] for dim in self.dims if dim in indexable_dims}
        for dim in indexes:
            if required_indexes[dim]: indexes[dim] = required_indexes[dim]

        # Generate Cartesian product
        keys = list(indexes.keys())
        values = list(indexes.values())

        for combo in product(*values):
            yield dict(zip(keys, combo))

    def load(self,indexes):
        if self.key is None: return
        ###
        ### Get data
        ###
        da = self.ds[self.key]
        self.data = da[indexes]
        self.indexes = indexes
        ###
        ### Work with a 2-d array da(lat,lon)
        ###
        if ["lat","lon"] != sorted(self.data.dims): 
            logger.info(f"Incorrect dimensions for {self.key}. Nothing to do")
            return

    def plot(self):
        if self.key is None: return

        fig, ax = self._create_map()
        self.ax = ax
        #
        autoScale = False
        #
        factor = self._get_factor()
        cmap   = self._get_colormap()
        levels = self._get_levels()
        ###
        ### Configure plot
        ###
        args_dict = {
                "cmap": cmap,
                "extend": 'max',
                "transform": crs.PlateCarree(),
                }
        if not autoScale: 
            args_dict['levels'] = levels
            args_dict['norm'] = BoundaryNorm(levels,cmap.N)
        ###
        ### Create plot
        ###
        fc = ax.contourf(self.data.lon,self.data.lat,factor*self.data,**args_dict)
        ###
        ### Generate colorbar
        ###
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(
            position='bottom',
            size="3%", 
            pad=0.4,
            axes_class=plt.Axes)
        cbar = fig.colorbar(fc, 
            orientation='horizontal',
            label = self._get_label(),
            cax=cax)
        cbar.set_ticks(levels[:-1])

    def add_title(self):
        ###
        ### Add title
        ###
        if self.key is None: return

        if 'time' in self.dims:
            it = self.indexes['time']
            time_fmt = self.ds.isel(time=it)['time'].dt.strftime(FMT_TIME).item()
            self.ax.set_title(f"valid: {time_fmt}", loc='right')
            self.ax.set_title(f"+{it:03d}h FCST", loc='left')

    def add_marker(self,lon,lat):
        ###
        ### Add vent location
        ###
        if self.key is None: return

        self.ax.plot(lon,lat,
            color='red',
            marker='^',
            zorder = 4,
            alpha = 0.5,
            transform=crs.PlateCarree())

    def add_logo(self,fname):
        ###
        ### Add logo
        ###
        if self.key is None: return

        image_data = Image.open(fname)
        image_box = OffsetImage(image_data, zoom=0.2)
        anno_box = AnnotationBbox(image_box, 
            xy=(0, 1), 
            xycoords='axes fraction', 
            box_alignment=(0, 1), 
            frameon=False)
        self.ax.add_artist(anno_box)
 
    def save(self):
        """
        Save a map as a png image
        """
        if self.key is None: return

        prefix = self.key
        suffix = ""

        if 'layer' in self.dims:
            ilayer = self.indexes['layer']
            suffix += f"_l{ilayer}"

        if 'time' in self.dims:
            it = self.indexes['time']
            suffix += f"_t{it:03}"

        fname = f"{prefix}{suffix}.png"
        logger.info(f"Saving {fname}")
        ##
        ## Save output file
        ##
        plt.savefig(fname,dpi=200,bbox_inches='tight')

    def _get_label(self):
        label = self.data.long_name.replace("_", " ")

        if label == 'tephra_grn_load': 
            label = "deposit thickness"

        if 'layer' in self.dims:
            i = self.indexes['layer']
            top = self.ds["layer_top"].isel(layer=i).item()
            bottom = self.ds["layer_bottom"].isel(layer=i).item()
            label = f"{label} (FL{bottom:.0f}-{top:.0f})"

        ### Append units
        units = self._get_units()
        if units:
            label = fr"{label} [${units}$]"
        else:
            label = f"{label}"
        return label

    def _get_units(self):
        if self.key in ["SO2_con_layer", "SO2_fl", "tephra_con_layer", "tephra_fl"]:
            units = 'mg~m^{-3}'
        elif self.key in ["tephra_col_mass", "tephra_col_mass_pm"]:
            units = 'g~m^{-2}'
        elif self.key in ["SO2_col_mass"]:
            units = 'DU'
        elif self.key in ["tephra_grn_load", "SO2_grn_load"]:
            units = 'mm'
        elif self.key in ["tephra_cloud_top", "SO2_cloud_top"]:
            units = 'km'
        else:
            units = '' 
        return units

    def _get_colormap(self):
        ###
        ### Define colormap
        ###
        if self.key in ["SO2_con_layer", "SO2_fl", "tephra_con_layer", "tephra_fl"]:
            cmap = ListedColormap(["cyan", "grey", "red"])
        elif self.key in ["tephra_col_mass", "tephra_grn_load", "SO2_grn_load"]:
            cmap = plt.cm.RdYlBu_r
        else:
            cmap = plt.cm.viridis
        return cmap

    def _get_levels(self):
        if self.key in ["SO2_con_layer", "SO2_fl", "tephra_con_layer", "tephra_fl"]:
            levels = [0.2, 2, 4, 6]
        elif self.key in ["tephra_col_mass"]:
            levels = [0.1, 0.2, 0.4, 1, 2, 4, 10, 20, 40, 100, 200]
        elif self.key in ["SO2_cloud_top", "tephra_cloud_top"]:
            levels = [4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
        else:
            levels = [0.1, 0.2, 0.4, 1, 2, 4, 10, 20]
        return levels

    def _get_factor(self):
        if self.key in ["SO2_con_layer", "SO2_fl", "tephra_con_layer", "tephra_fl"]:
            factor = 1E3
        elif self.key in ["SO2_cloud_top", "tephra_cloud_top"]:
            factor = 1E-3
        else:
            factor = 1.0
        return factor
 
    def _create_map(self):
        proj = crs.PlateCarree()
        fig, ax = plt.subplots( subplot_kw={'projection': proj} )
        ###
        ### Add map features
        ###
        BORDERS = cfeature.NaturalEarthFeature(
                scale     = '10m' if self.hires else '50m',
                category  = 'cultural',
                name      = 'admin_0_countries',
                edgecolor = 'gray',
                facecolor = 'none'
                )
        LAND = cfeature.NaturalEarthFeature(
                category  = 'physical',
                name      = 'land',
                scale     = '10m' if self.hires else '50m',
                edgecolor = 'none',
                facecolor = 'lightgrey',
                alpha     = 0.8
                )
        ax.add_feature(LAND, zorder=0)
        ax.add_feature(BORDERS, linewidth=0.4)
        ###
        ### Add grid lines
        ###
        gl = ax.gridlines(
            crs         = crs.PlateCarree(),
            draw_labels = True,
            linewidth   = 0.5,
            color       = 'gray',
            alpha       = 0.5,
            linestyle   = '--')
        gl.top_labels    = False
        gl.right_labels  = False
        gl.ylabel_style  = {'rotation': 90}
        return (fig,ax)

def main(args):
    addVolcano = True
    addLogo    = True
    #
    if args.lat is None or args.lon is None: addVolcano = False
    if args.logo is None: addLogo = False
    #
    c = MapImage(args.netcdf)
    c.set_key(args.key)
    for ind in c.key_indexes(time = args.times, layer = args.layers):
        c.load(ind)
        c.plot()
        if addVolcano: c.add_marker(args.lon,args.lat)
        if addLogo: c.add_logo(args.logo)
        c.add_title()
        c.save()

if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Plot a map from a NetCDF file")
    parser.add_argument("--netcdf", metavar='file',      type=str,   help="Path to NetCDF file", required=True)
    parser.add_argument("--key",    metavar='variable',  type=str,   help="Variable name in NetCDF", required=True)
    parser.add_argument('--lat',    metavar='latitude',  type=float, help='Volcano latitude')
    parser.add_argument('--lon',    metavar='longitude', type=float, help='Volcano longitude')
    parser.add_argument("--times",  metavar='indexes',   type=int,   help='List of time step index', nargs='+')
    parser.add_argument("--layers", metavar='indexes',   type=int,   help='List of layer indexes', nargs='+')
    parser.add_argument("--logo",   metavar='file',      type=str,   help='Logo image')
    args = parser.parse_args()

    logger.info("Creating map output file...")
    main(args)
    logger.info("Done!")
