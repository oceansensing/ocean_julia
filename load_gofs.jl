using PyCall, NaNMath, Dates, Missings, Glob, Interpolations
using GibbsSeaWater, NetCDF, NCDatasets

gofs_url = "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0";
gofs = NCDataset(gofs_url, "r");

latmin, latmax = 30.0, 43.0;
lonmin, lonmax = -82.0, -45.0;

tnow = Dates.datetime2unix(Dates.now());
t0 = tnow;

lat_gofs = gofs["lat"][:];
lon_gofs = gofs["lon"][:];
t_gofs = Dates.datetime2unix.(gofs["time"][:]);

(t0min, t0ind) = findmin(abs.(tnow .- t_gofs));
