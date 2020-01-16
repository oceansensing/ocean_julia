# Install Julia and enter the directory where this script is located:
# cd("/Users/Shabangin/Research/NPP/")
#
# run this script by callling:
# include("npp_pkg_setup.jl")
#
# Note: You should only need to run this script ONCE for each Julia installation, although it doesn't hurt anything to run it more than once. The first time you run it it'll take a long while.
#
# 2019-11-19: Donglai Gong
# 2019-12-14: Donglai Gong, clean it up so that packages not needed are not installed.

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

using Pkg
#Pkg.add("Debugger")
Pkg.add("IJulia")
Pkg.add("Conda")
Pkg.add("PyCall")
Pkg.add("Dates")
Pkg.add("Statistics")
Pkg.add("NaNMath")
#Pkg.add("RandomForests")
Pkg.add("DataFrames")
Pkg.add("Missings")
Pkg.add("NCDatasets")
Pkg.add(PackageSpec(url="https://github.com/kouketsu/GSW.jl", rev="master"))
Pkg.add(PackageSpec(url="https://github.com/JuliaGeo/GeoDatasets.jl", rev="master"))

Pkg.add("JLD2")
Pkg.add("FileIO")
Pkg.add("Glob")
Pkg.add("CSV")

Pkg.add("PyPlot")
Pkg.add("Pandas")
Pkg.add("WebIO")

using Conda
Conda.add("cmocean",channel="conda-forge")

if Sys.iswindows() == false
    ENV["PYTHON"] = ENV["HOME"] * "/anaconda3/bin/python3"
    #ENV["PYTHON"] = ENV["HOME"] * "/.julia/conda/3/bin/python3"
else
    ENV["PYTHON"] = "C:" * ENV["HOMEPATH"] * "\\.julia\\conda\\3\\python.exe"
end
Pkg.build("PyCall")

using CSV, Glob, FileIO, JLD2, DataFrames #Serialization, FileIO, JLD
using Missings, Dates
using GSW, Statistics, NaNMath
using IJulia, PyPlot, PyCall

Pkg.update()
