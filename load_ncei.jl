# setting up needed Julia packages
using Pkg
Pkg.add(["HTTP","JSON2"])
using HTTP, JSON2

# https://www.ncdc.noaa.gov/data-access
# define NOAA NDBC data access web address & data access token
base_ncei_url = "https://www.ncei.noaa.gov/access/services/data/v1";
#token = "FMAPxDJFXSsMKsJWUDiOhZkINCfYGZqX"; # need to request your own from https://www.ncdc.noaa.gov/cdo-web/token

dataset = "data";
datasetid = "GHCND"
datatypeid = "TOBS"
stationid = Dict(["williamsburg" => "GHCND:USC00449151", "newportnews" => "GHCND:USW00093741"]) # Williamsburg, Newport News
startdate = "2020-04-01"
enddate = "2020-04-03"
rd = HTTP.get(base_ncei_url * endpoint * "?dataset=" * datasetid * "&stationid=" * stationid["williamsburg"] * "&datatypeid=" * datatypeid * "&startdate=" * startdate * "&enddate=" * enddate, ["token" => token])
@pretty String(rd.body)
