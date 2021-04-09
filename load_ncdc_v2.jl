# setting up needed Julia packages
using Pkg
Pkg.add(["HTTP","JSON2"])
using HTTP, JSON2

# define NOAA NDBC data access web address & data access token
base_ncdc_url = "https://www.ncdc.noaa.gov/cdo-web/api/v2/";
token = "FMAPxDJFXSsMKsJWUDiOhZkINCfYGZqX"; # need to request your own from https://www.ncdc.noaa.gov/cdo-web/token

# use HTTP.get to obtain location categories, see https://www.ncdc.noaa.gov/cdo-web/webservices/v2#gettingStarted
endpoint = "locationcategories";
r = HTTP.get(base_ncdc_url * endpoint, ["token" => token])
@pretty String(r.body)

endpoint = "datasets";
id = "/GHCND"
rds = HTTP.get(base_ncdc_url * endpoint * id, ["token" => token])
@pretty String(rds.body)

endpoint = "data";
datasetid = "GHCND"
datatypeid = "TOBS"
stationid = Dict(["williamsburg" => "GHCND:USC00449151", "newportnews" => "GHCND:USW00093741"]) # Williamsburg, Newport News
startdate = "2020-04-01"
enddate = "2020-04-03"
rd = HTTP.get(base_ncdc_url * endpoint * "?datasetid=" * datasetid * "&stationid=" * stationid["williamsburg"] * "&datatypeid=" * datatypeid * "&startdate=" * startdate * "&enddate=" * enddate, ["token" => token])
@pretty String(rd.body)
