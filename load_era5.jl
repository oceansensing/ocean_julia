using Pkg

using PyCall, Glob, NCDatasets, Missings, Interpolations, NaNMath, Dates, Zarr, Plots
using C2PO

cmo = pyimport("cmocean")
gsw = pyimport("gsw")
netCDF4 = pyimport("netCDF4")

mutable struct ERA5surf
    time::Array{AbstractFloat}
    longitude::Array{AbstractFloat}
    latitude::Array{AbstractFloat}
    u10::Array{AbstractFloat, 3} # u wind at 10 m
    v10::Array{AbstractFloat, 3} # v wind at 10 m
    t2m::Array{AbstractFloat, 3} # temp at 2 m
    d2m::Array{AbstractFloat, 3} # dew point at 2 m
    msl::Array{AbstractFloat, 3} # mean sea level
    mwd::Array{AbstractFloat, 3} # mean wave direction
    mwp::Array{AbstractFloat, 3} # mean wave period
    sst::Array{AbstractFloat, 3} # sea surface temperature
    swh::Array{AbstractFloat, 3} # significant wave height
    sp::Array{AbstractFloat, 3} # surface air pressure
    tp::Array{AbstractFloat, 3} # total precipation
end

plotly()

era5dir = "/Users/gong/oceansensing Dropbox/C2PO/Data/ERA5/"
yyyy = 2019;
mm = 9;

datadir = era5dir * string(yyyy) * "/";
era5file = "ERA5_surface_hourly_" * string(yyyy) * string(mm, pad = 2) * ".nc"

era5 = NCDataset(datadir * era5file, "r");

time = era5["time"][:];
lon = era5["longitude"][:];
lat = era5["latitude"][:];