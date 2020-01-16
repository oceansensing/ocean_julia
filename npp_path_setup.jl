# this script contains the directory and filenames that other scripts depends, please update if necessary
# 2019-12-12, Donglai Gong: I revamped this script so that it can automatically determine whether it's Donglai's, Mirella's, or Tristan's computer and setup the directories for each accordingly. Please let me know if you encounter issue with this script.
using Glob

if Sys.iswindows() == false
    dirstr = split(pwd(),"/");
    if "gong" in dirstr
        println("Donglai's system.")
        ctddir = string(ENV["HOME"], "/Research/NPP/");
        ctddatadir = string(ctddir,"data/CTD/cnv3/");
        ctddatapath = Glob.glob("*.cnv",ctddatadir);

        ctddatadircsv = string(ctddir,"data/CTD/csv3/");
        ctddatapathcsv = Glob.glob("*.csv",ctddatadircsv);

        workdir = string( ctddir, "data/CTD/");
        figoutdir = string(ctddir, "figures3/");
        filepath = string(workdir, "NPP_ctd_data", ".jld2");
        ctdlocationpath = string(workdir, "NPP_ctd_locations.csv");
    elseif "Shabangin" in dirstr
        println("Mirella's system.")
        ctddir = string(ENV["HOME"], "/Research/NPP/"); # this is where your data directory is specified

        ctddatadir = string(ctddir,"cnv3/"); # this is where your CNV data directory is specified
        ctddatapath = Glob.glob("*.cnv",ctddatadir);

        ctddatadircsv = string(ctddir,"csv3/") # this is your your CSV data directory is specified (not needed if just using CNV files)
        ctddatapathcsv = Glob.glob("*.csv",ctddatadircsv);

        workdir = ctddir; # this is where your work directory is specified
        filepath = string(workdir, "NPP_ctd_data", ".jld2"); # define the name of the data file to store the loaded CTD data in JLD2 format
        figoutdir = string(ctddir, "figures3/");
        ctdlocationpath = string(workdir, "NPP_ctd_locations.csv");
    end
else
    println("Tristan's system.")
    ctddir = string("C:", ENV["HOMEPATH"], "\\Documents\\Research\\NPP\\"); # this is where your data directory is specified

    ctddatadir = string(ctddir,"cnv3\\"); # this is where your CNV data directory is specified
    ctddatapath = Glob.glob("*.cnv",ctddatadir);

    ctddatadircsv = string(ctddir,"csv3\\"); # this is your your CSV data directory is specified (not needed if just using CNV files)
    ctddatapathcsv = Glob.glob("*.csv",ctddatadircsv);

    workdir = ctddir; # this is where your work directory is specified
    filepath = string(workdir, "NPP_ctd_data", ".jld2"); # define the name of the data file to store the loaded CTD data in JLD2 format
    figoutdir = string(ctddir, "figures3\\");
    ctdlocationpath = string(workdir, "NPP_ctd_locations.csv");
end
