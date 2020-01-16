#ENV["JULIA_PARDISO"] = "/Users/gong/lib/Pardiso/libpardiso600-MACOS-X86-64.dylib"
ENV["PYTHON"] = string(ENV["HOME"], "/anaconda3/bin/python3")

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

using Pkg
