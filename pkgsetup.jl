# 2020-02-28: Donglai Gong

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

using Pkg
Pkg.add("IJulia")
Pkg.add("ScikitLearn")
Pkg.add("Measurements")
Pkg.add("MLJ")
Pkg.add("Flux")
Pkg.add("LsqFit")
Pkg.add("TimeSeries")
Pkg.add("Indicators")
Pkg.add("DynamicalSystems")
Pkg.add("Oceananigans")
Pkg.add("DifferentialEquations")
Pkg.add("Quaternions")
Pkg.add("GLM")
Pkg.add("Convex")
Pkg.add("Glob")
Pkg.add("GeoJSON")
Pkg.add("Pandas")
Pkg.add("ONNX")
Pkg.add("JuliaDB")
Pkg.add("MATLAB")
Pkg.add("TextAnalysis")
Pkg.add("Makie")
Pkg.add("StatsPlots")
Pkg.add("Interact")
Pkg.add(PackageSpec(url="https://github.com/kouketsu/GSW.jl", rev="master"))

using Conda
Conda.add("cmocean",channel="conda-forge")

if Sys.iswindows() == false
    #ENV["PYTHON"] = ENV["HOME"] * "/anaconda3/bin/python3"
    ENV["PYTHON"] = ENV["HOME"] * "/.julia/conda/3/bin/python3"
else
    ENV["PYTHON"] = "C:" * ENV["HOMEPATH"] * "\\.julia\\conda\\3\\python.exe"
end
Pkg.build("PyCall")

using NCDatasets
using CSV, Glob, FileIO, JLD2, DataFrames #Serialization, FileIO, JLD
using Missings, Dates
using GSW, Statistics, NaNMath
using IJulia, PyPlot, PyCall
using Makie, Interact
