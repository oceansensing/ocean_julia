ENV["JULIA_PARDISO"] = "/Users/gong/lib/Pardiso/libpardiso600-MACOS-X86-64.dylib"
ENV["PYTHON"] = string(ENV["HOME"], "/anaconda3/bin/python3")
#ENV["PYTHON"] = string(ENV["HOME"], "/.julia/conda/3/bin/python3")
#ENV["CONDA_JL_HOME"] =

ocean_julia = ENV["HOME"] * "/GitHub/ocean_julia"
if (ocean_julia in LOAD_PATH) == false
    push!(LOAD_PATH, ocean_julia);
end

ocean_julia_npp = ENV["HOME"] * "/GitHub/ocean_julia/NPP"
if (ocean_julia_npp in LOAD_PATH) == false
    push!(LOAD_PATH, ocean_julia_npp);
end

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

using Pkg
