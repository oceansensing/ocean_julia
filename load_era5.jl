using Pkg

using PyCall, Glob, NCDatasets, Missings, Interpolations, NaNMath, Dates, Zarr
using Plots

cmo = pyimport("cmocean")
gsw = pyimport("gsw")
netCDF4 = pyimport("netCDF4")
plotly()

era5dir = ""