using PyCall, NaNMath, Dates, Missings, Glob, Interpolations
using GibbsSeaWater, NetCDF, NCDatasets

gofs_url = "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0";
gofs = NCDataset(gofs_url, "r");

lat_gofs = gofs["lat"][:];
lon_gofs = gofs["lon"][:];
t_gofs = gofs["time"][:];

