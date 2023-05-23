using Pkg, PyCall, NaNMath, Dates, Missings, Glob, Interpolations
using GibbsSeaWater, NetCDF, NCDatasets
import C2PO



# Open HYCOM GOFS 3.1's NetCDF OpenDAP data link
gofs_url = "https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0";
gofs = NCDataset(gofs_url, "r");

# define study region's box and time of interest
latmin, latmax = 30.0, 43.0;
lonmin, lonmax = -82.0, -45.0;
tnow = Dates.datetime2unix(Dates.now());
t0 = tnow;

# extract model grid and time indices
lat_gofs = gofs["lat"][:];
lon_gofs = gofs["lon"][:] - 360.0;
t_gofs = Dates.datetime2unix.(gofs["time"][:]);

# find the grid and time indices of interesting
(t0min, t0ind) = findmin(abs.(tnow .- t_gofs));
latind = findall(latmin .<= lat_gofs .<= latmax);
lonind = findall(lonmin .<= lon_gofs .<= lonmax);
