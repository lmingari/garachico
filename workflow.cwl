#!/usr/bin/env cwl-runner
cwlVersion: v1.2
$namespaces:
  s: https://schema.org/
schemas:
- http://schema.org/version/9.0/schemaorg-current-http.rdf

s:softwareVersion: 0.0.1

s:author:
  - class: s:Person
    s:name: Leonardo Mingari
    s:email: mailto:lmingari@csic.es
    s:identifier: https://orcid.org/0000-0002-6584-4699

s:codeRepository: https://https://github.com/lmingari/garachico

$graph:
  - id: config
    class: CommandLineTool
    baseCommand: ["fill_template.py"]
    inputs:
      template:
        type: string
        inputBinding: {prefix: --template}
      initial_condition:
        type: "#main/InitialCondition"
        inputBinding: {prefix: --INITIAL}
      meteo_database:
        type: "#main/MeteoDatabase"
        inputBinding: {prefix: --METEO_DATABASE}
      meteo:
        type: string
        inputBinding: {prefix: --METEO_FILE}
      dictionary:
        type: string?
        inputBinding: {prefix: --METEO_DICTIONARY}
      restart:
        type: string?
        inputBinding: {prefix: --RESTART_FILE}
      levels:
        type: string?
        inputBinding: {prefix: --LEVELS_FILE}
      date:
        type: string?
        inputBinding: {prefix: --date}
      start_time:
        type: float
        inputBinding: {prefix: --START_TIME}
      end_time:
        type: float
        inputBinding: {prefix: --END_TIME}
      west_lon:
        type: float
        inputBinding: {prefix: --LONMIN}
      east_lon:
        type: float
        inputBinding: {prefix: --LONMAX}
      north_lat:
        type: float
        inputBinding: {prefix: --LATMAX}
      south_lat:
        type: float
        inputBinding: {prefix: --LATMIN}
      dx:
        type: float
        inputBinding: {prefix: --DX}
      dy:
        type: float
        inputBinding: {prefix: --DY}
      nlevels:
        type: int
        inputBinding: {prefix: --NZ}
      lon_vent:
        type: float
        inputBinding: {prefix: --LONV}
      lat_vent:
        type: float
        inputBinding: {prefix: --LATV}
      elevation_vent:
        type: float
        inputBinding: {prefix: --ELEVATION}
      particle_cutoff:
        type: float?
        inputBinding: {prefix: --PARTICLE_CUTOFF}
    outputs:
      inp:
        type: File
        outputBinding:
          glob: "*.inp"

  - id: gfs_times
    class: CommandLineTool
    baseCommand: ["gfs_times.py"]
    inputs:
      reference_date:
        type: string
        inputBinding: {prefix: --date}
      hmin:
        type: float
        inputBinding: {prefix: --hinit}
      hmax:
        type: float
        inputBinding: {prefix: --hend}
      step:
        type: int?
        inputBinding: {prefix: --step}
    outputs:
      cycle:
        type: int
        outputBinding:
          glob: gfs.json
          loadContents: true
          outputEval: $(JSON.parse(self[0].contents).cycle)
      date:
        type: string
        outputBinding:
          glob: gfs.json
          loadContents: true
          outputEval: $(JSON.parse(self[0].contents).date)
      tmin:
        type: int
        outputBinding:
          glob: gfs.json
          loadContents: true
          outputEval: $(JSON.parse(self[0].contents).tmin)
      tmax:
        type: int
        outputBinding:
          glob: gfs.json
          loadContents: true
          outputEval: $(JSON.parse(self[0].contents).tmax)
    stdout: gfs.json
    requirements:
      InlineJavascriptRequirement: {}

  - id: set_emission
    class: CommandLineTool
    label: Generate eruptive phases file for FALL3D
    baseCommand: echo
    arguments: ["Creating phases file"]
    inputs:
      elevation_vent: float
      phases:
        type:
          type: array
          items: "#main/PhaseDefinition"
    outputs: 
      table:
        label: Emission file for FALL3D
        type: File
        outputBinding:
          glob: "phases.dat"
    requirements:
      InitialWorkDirRequirement:
        listing:
          - entryname: "phases.dat"
            entry: |
              ${
                  let height_avl;
                  let lines = [];
                  for (let p of inputs.phases) {
                      height_avl = (p.height > inputs.elevation_vent) ? (p.height - inputs.elevation_vent) : 0;
                      lines.push(p.start.toFixed(2) + "\t" + 
                                 p.end.toFixed(2)   + "\t" + 
                                 height_avl.toFixed(2) );
                  }
                  return lines.join("\n");
              }

  - id: gfs
    class: CommandLineTool
    baseCommand: ["gfs.py"]
    inputs:
      date:
        type: string
        inputBinding: {prefix: --date}
      lonmin:
        type: float
        inputBinding:
          prefix: --lon
          position: 1
      lonmax:
        type: float
        inputBinding: 
          position: 2
      latmin:
        type: float
        inputBinding: 
          prefix: --lat
          position: 3
      latmax:
        type: float
        inputBinding: 
          position: 4
      tmin:
        type: int
        inputBinding: 
          prefix: --time
          position: 5
      tmax:
        type: int
        inputBinding: 
          position: 6
      resolution:
        type: float
        inputBinding: {prefix: --res}
      cycle:
        type: int
        inputBinding: {prefix: --cycle}
      step:
        type: int
        inputBinding: {prefix: --step}
    outputs:
      gribs:
        type: File[]
        outputBinding:
          glob: "*.f???"

  - id: gfs_table
    class: ExpressionTool
    inputs:
      resolution: float
    outputs:
      nlev: int
      levels: float[]
    expression: |
      ${
          var selectedList;
          var res = inputs.resolution;
          if (res === 0.5) {
              selectedList = [1000,975,950,925,900,875,850,825,800,775,750,725,700,675,650,625,600,575,550,525,500,475,450,425,400,375,350,325,300,275,250,225,200,175,150,125,100,70,50,40,30,20,15,10,7,5,3,2,1,0.7,0.4,0.2,0.1,0.07,0.04,0.02,0.01];
          } else {
              selectedList = [1000,975,950,925,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100,70,50,40,30,20,15,10,7,5,3,2,1,0.7,0.4,0.2,0.1,0.07,0.04,0.02,0.01];
          }
          return { 
              levels: selectedList,
              nlev: selectedList.length
          };
      }
    requirements:
      InlineJavascriptRequirement: {}

  - id: grib2nc
    class: CommandLineTool
    baseCommand: sh
    arguments: ["grib2nc.sh"]
    inputs:
      gribs:
        type: File[]
        inputBinding: {position: 1}
      nlev: int
      levels: float[]
    outputs:
      netcdf:
        type: File
        outputBinding:
          glob: "*.nc"
    stdout: stdout.wgrib
    requirements:
      InlineJavascriptRequirement: {}
      InitialWorkDirRequirement:
        listing:
          - entryname: table.gfs
            entry: |
              $lev_type  100:lev:pressure level:mb:0.01
              $nlev      $(inputs.nlev)
              $levs      $(inputs.levels.join())
          - entryname: grib2nc.sh
            entry: |
              #!/bin/bash
              variables="HGT|TMP|RH|UGRD|VGRD|VVEL|PRES|PRATE|LAND|HPBL|SFCR"
              levels="([0-9]*[.])?[0-9]+ mb|surface|2 m above ground|10 m above ground"
              for file in "$@"; do
                  echo "Processing $file"
                  wgrib2 $file \
                    -match ":(\${variables}):" \
                    -match ":(\${levels}):" \
                    -nc_table "table.gfs" \
                    -append \
                    -nc3 \
                    -netcdf \
                    gfs.nc                    
              done

  - id: runner
    class: CommandLineTool
    label: Run FALL3D model
    baseCommand: mpirun
    arguments:
      - prefix: -n
        valueFrom: $(inputs.nx * inputs.ny * inputs.nz)
      - "Fall3d.x"
    inputs:
      task:
        type: string
        default: all
        inputBinding: {position: 0}
      inp: 
        type: File
        inputBinding: {position: 1}
      nx:
        type: int
        inputBinding: {position: 2}
      ny:
        type: int
        inputBinding: {position: 3}
      nz:
        type: int
        inputBinding: {position: 4}
      meteo: File
      emission: File
    outputs:
      stdout:
        type: stdout
      stderr:
        type: stderr
      log:
        type: File
        outputBinding:
          glob: "*.Fall3d.log"
      res:
        type: File
        outputBinding:
          glob: "*.res.nc"
    stdout: fall3d.out
    stderr: fall3d.err
    requirements:
      InlineJavascriptRequirement: {}
      InitialWorkDirRequirement:
        listing:
          - $(inputs.inp)
          - $(inputs.meteo)
          - $(inputs.emission)

  - id: plotter
    class: CommandLineTool
    baseCommand: ["plotter.py"]
    inputs:
      key:
        type: string
        inputBinding: {prefix: --key}
      netcdf:
        type: File
        inputBinding: {prefix: --netcdf}
      lat:
        type: float?
        inputBinding: {prefix: --lat}
      lon:
        type: float?
        inputBinding: {prefix: --lon}
      times:
        type: int[]?
        inputBinding: {prefix: --times}
      layers:
        type: int[]?
        inputBinding: {prefix: --layers}
      logo:
        type: string?
        inputBinding: {prefix: --logo}
    outputs:
      pngs:
        type: File[]
        outputBinding: {glob: "*.png"}   

  - id: cog_creator
    class: CommandLineTool
    label: Generate a tif output from a netcdf
    baseCommand: 'createCOG.py'
    arguments: []
    doc: >
      This tool reads a variable (key) from a FALL3D output 
      file in netCDF format (netcdf) and produces a single
      Cloud-Optimized GeoTIFF (COG) file in geotiff format
    inputs:
      netcdf:
        label: FALL3D output file in netCDF format
        type: File
        inputBinding: {prefix: --fname}
      key:
        label: Variable key name to be plotted
        type: string
        inputBinding: {prefix: --key}
      time:
        label: Time step (zero-based indexed)
        type: int
        inputBinding: {prefix: --time}
    outputs:
      cog:
        label: Georeferenced image in GeoTIFF file format
        type: File
        outputBinding:
          glob: "*.tif"

  - id: gather_files
    class: ExpressionTool
    inputs:
      folder: string
      files: File[]
    outputs:
      collected_folder: Directory
    expression: |
      ${
          var dir = {
              "class": "Directory",
              "basename": inputs.folder,
              "listing": inputs.files
          };
          return {"collected_folder": dir};
      }
    requirements:
      InlineJavascriptRequirement: {}

  - id: main
    class: Workflow
    label: FALL3D forecast driven by GFS
    inputs:
      volcano: string
      lat_vent: float
      lon_vent: float
      elevation_vent: float
      date: string
      start_time: float
      end_time: float
      particle_cutoff: float
      meteo_database: 
        type: "#main/MeteoDatabase"
      initial_condition:
        type: "#main/InitialCondition"
      restart: string?
      west_lon: float
      east_lon: float
      north_lat: float
      south_lat: float
      dx: float
      dy: float
      nlevels: int
      # Eruptive phases
      phases:
        type:
          type: array
          items: "#main/PhaseDefinition"
      # GFS inputs
      resolution: float
      step: int
      # MPI inputs
      nx: int
      ny: int
      nz: int
      # Plot inputs
      times: int[]
      keys: string[]
      layers: int[]?
    outputs: 
      pngs:
        type: Directory
        outputSource: png_folder/collected_folder
      cogs:
        type: Directory
        outputSource: cog_folder/collected_folder
    steps:
      configure:
        run: "#config"
        in:
          template: {default: "/app/template.inp"}
          meteo: 
            source: gfs_netcdf/netcdf
            valueFrom: $(self.basename)
          dictionary: {default: "/app/GFS.tbl"}
          meteo_database: meteo_database
          date: date
          start_time: start_time
          end_time: end_time
          west_lon: west_lon
          east_lon: east_lon
          north_lat: north_lat
          south_lat: south_lat
          dx: dx
          dy: dy
          nlevels: nlevels
          initial_condition: initial_condition
          restart: restart
          lon_vent: lon_vent
          lat_vent: lat_vent
          elevation_vent: elevation_vent
          particle_cutoff: particle_cutoff
        out: [inp]
      emission:
        run: "#set_emission"
        in:
          elevation_vent: elevation_vent
          phases: phases
        out: [table]
      gfs_times:
        run: "#gfs_times"
        in:
          reference_date: date
          hmin: start_time
          hmax: end_time
          step: step
        out: [cycle,date,tmin,tmax]
      gfs_gribs:
        run: "#gfs"
        in:
          date: gfs_times/date
          tmin: gfs_times/tmin
          tmax: gfs_times/tmax
          cycle: gfs_times/cycle
          lonmin: 
            source: west_lon
            valueFrom: $(self - inputs.resolution)
          lonmax: 
            source: east_lon
            valueFrom: $(self + inputs.resolution)
          latmin:
            source: south_lat
            valueFrom: $(self - inputs.resolution)
          latmax:
            source: north_lat 
            valueFrom: $(self + inputs.resolution)
          resolution: resolution
          step: step
        out: [gribs]
      gfs_levels:
        run: "#gfs_table"
        in:
          resolution: resolution
        out: [nlev,levels]
      gfs_netcdf:
        run: "#grib2nc"
        in:
          gribs: gfs_gribs/gribs
          nlev: gfs_levels/nlev
          levels: gfs_levels/levels
        out: [netcdf]
      run_fall3d:
        run: "#runner"
        in:
          inp: configure/inp
          nx: nx
          ny: ny
          nz: nz
          meteo: gfs_netcdf/netcdf
          emission: emission/table
        out: [log,res]
      plot_maps:
        run: "#plotter"
        scatter: key
        in:
          key: keys
          times: times
          layers: layers
          netcdf: run_fall3d/res
          lon: lon_vent
          lat: lat_vent
          logo: {default: "/app/logo.png"}
        out: [pngs]
      create_cogs:
        run: "#cog_creator"
        scatter: [time,key]
        scatterMethod: flat_crossproduct
        in:
          netcdf: run_fall3d/res
          key: keys
          time: times
        out: [cog]
      png_folder:
        run: "#gather_files"
        in:
          files: 
            source: plot_maps/pngs
            valueFrom: $(self.flat())
          folder:
            source: [date,volcano]
            valueFrom: $(self.join("/"))
        out: [collected_folder]
      cog_folder:
        run: "#gather_files"
        in:
          files: create_cogs/cog
          folder:
            source: [date,volcano]
            valueFrom: $(self.join("/"))
        out: [collected_folder]
    requirements:
      StepInputExpressionRequirement: {}
      InlineJavascriptRequirement: {}
      ScatterFeatureRequirement: {}
      MultipleInputFeatureRequirement: {}
      NetworkAccess:
        networkAccess: true
      DockerRequirement:
        dockerPull: docker.io/lmingari/garachico:1.2
      SchemaDefRequirement:
        types:
          - name: MeteoDatabase
            type: enum
            symbols:
              - GFS
              - WRF
              - ERA5
          - name: InitialCondition
            type: enum
            symbols:
              - NONE
              - RESTART
              - INSERTION
          - name: PhaseDefinition
            type: record
            fields:
              start: float
              end: float
              height: float
