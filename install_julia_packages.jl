Pkg.activate(".julia/environments/ocean")
Pkg.add("DataFrames")
Pkg.add("BenchmarkTools")
Pkg.add("NaNMath")
Pkg.add("Missings")
Pkg.add("Glob")
Pkg.add("CSV")
Pkg.add("Pluto")
Pkg.add("DifferentialEquations")
Pkg.add("Optimization")
Pkg.add("NCDatasets") # installs  [eb928a42] + prrte_jll v3.0.2+0
Pkg.add(url="https://github.com/TEOS-10/GibbsSeaWater.jl", rev="master")
Pkg.add("AppleAccelerate")
Pkg.add("PyCall")
#Pkg.add("Flux") # this is causing precompiling issues


#Pkg.add(url="https://github.com/Alexander-Barth/ROMS.jl", rev="master") # this causes precompile issues, not working with Julia 1.10