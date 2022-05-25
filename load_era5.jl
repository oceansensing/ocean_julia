using Pkg

using PyCall, Glob, NCDatasets, Missings, Interpolations, NaNMath, Dates, Zarr
using Plots, PlotlySave

cmo = pyimport("cmocean")
gsw = pyimport("gsw")
netCDF4 = pyimport("netCDF4")
plotly()

era5dir = "/Users/gong/oceansensing Dropbox/C2PO/Data/ERA5/"
yyyy = 2019;

datadir = era5dir * string(yyyy) * "/";