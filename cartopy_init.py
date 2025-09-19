import cartopy.feature as cfeature

for resolution in ["10m","50m"]:
    land = cfeature.NaturalEarthFeature(
            "physical", 
            "land", 
            resolution)
    borders = cfeature.NaturalEarthFeature(
            "cultural", 
            "admin_0_countries", 
            resolution)

    _ = list(land.geometries())
    _ = list(borders.geometries()) 
