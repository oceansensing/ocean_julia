ENV["PYTHON"] = "/Users/gong/opt/anaconda3/bin/python"

seaexplorer = ENV["HOME"] * "/GitHub/jlglider/seaexplorer"
if (seaexplorer in LOAD_PATH) == false
    push!(LOAD_PATH, seaexplorer);
end

slocum = ENV["HOME"] * "/GitHub/jlglider/slocum"
if (slocum in LOAD_PATH) == false
    push!(LOAD_PATH, slocum);
end

#if (pwd() in LOAD_PATH) == false
#    push!(LOAD_PATH, pwd());
#end
