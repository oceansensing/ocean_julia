#try
#    using Revise
#catch e
#    @warn "Error initializing Revise" exception=(e, catch_backtrace())
#end
ENV["PYTHON"] = "/Users/gong/anaconda3/bin/python"

ocean_julia_path = ENV["HOME"] * "/GitHub/jlglider/ocean_julia"
if (ocean_julia_path in LOAD_PATH) == false
    push!(LOAD_PATH, ocean_julia_path);
end

seaexplorer_path = ENV["HOME"] * "/GitHub/jlglider/seaexplorer"
if (seaexplorer_path in LOAD_PATH) == false
    push!(LOAD_PATH, seaexplorer_path);
end

slocum_path = ENV["HOME"] * "/GitHub/jlglider/slocum"
if (slocum_path in LOAD_PATH) == false
    push!(LOAD_PATH, slocum_path);
end

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end
