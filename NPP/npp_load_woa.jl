module WOA

using NCDatasets, Glob, CSV, Missings
using PyPlot
plt = PyPlot;

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

#import C2PO: gc_distance, rad2deg, deg2rad, histc, meshgrid, nan, findNaNmin, findNaNmax # this module stores additional function needed for analysis (e.g. gc_distance)

if (@isdefined ctd) == false
    import npp_load_ctd: ctd, npp
end

woadir = string(ENV["HOME"], "/Research/data/WOA18/");
woatemppath = Glob.glob("woa*81B0_t*.nc", woadir);
woasaltpath = Glob.glob("woa*81B0_s*.nc", woadir);
ctdlocationpath = string(ENV["HOME"], "/Research/NPP/data/CTD/NPP_ctd_locations.csv")
figoutdir = string(ENV["HOME"], "/Research/NPP/figures_WOA_NPP/")

# location NPP CTD location data
ctdloc = CSV.read(ctdlocationpath); # reading in the NPP CTD location data
lonctd = convert(Vector,ctdloc.lon); # converting from dataframe to vector
latctd = convert(Vector,ctdloc.lat);

mutable struct Profile
    z::Array{Float64}
    temp::Array{Float64}
    salt::Array{Float64}
end

mutable struct Season
    name::String
end

season = Dict([(1,"Winter"), (2,"Spring"), (3,"Summer"), (4,"Fall")])

# load WOA18 location data for each season
for i = 1:length(woatemppath)
    println(season[i])
    dstemp = Dataset(woatemppath[i]);
    dssalt = Dataset(woasaltpath[i]);
    latnc = dstemp["lat"];
    lonnc = dstemp["lon"];
    depthnc = dstemp["depth"];
    lon = Float64.(collect(Missings.replace(lonnc[:],NaN)));
    lat = Float64.(collect(Missings.replace(latnc[:],NaN)));
    z = -Float64.(collect(Missings.replace(depthnc[:],NaN)));
    #(LON,LAT) = meshgrid(lon,lat);

    prof = Profile[];

    tempnc = dstemp["t_an"];
    saltnc = dssalt["s_an"];
    temp = Float64.(collect(Missings.replace(tempnc[:],NaN)));
    salt = Float64.(collect(Missings.replace(saltnc[:],NaN)));

    for j = 1:length(lonctd)
        latind = findmin(abs.(lat .- latctd[j]))[2];
        lonind = findmin(abs.(lon .- lonctd[j]))[2];;
        push!(prof,Profile(z,temp[lonind,latind,:,1],salt[lonind,latind,:,1]));
    end

    for j = 1:length(ctd)
        println(string("Sta ", j))
        gind = findall(ctd[j].dataflag .> 0);

        figT = plt.figure(1)
        clf()
        axT = figT.add_subplot(1, 1, 1)

        axT.plot(prof[j].temp, prof[j].z, c="black")
        axT.plot(ctd[j].temp[gind], ctd[j].z[gind], c="red")
        axT.set_title(string("Station ", j, " (Temperature)"))
        axT.set_xlabel("Temperature (C)")
        axT.set_ylabel("Depth (m)")
        axT.legend(["WOA", "NPP"])
        plt.savefig(string(figoutdir, "temp_sta_", string(j/10000)[end-1:end], "_", season[i], ".png"), dpi = 300)

        figS = plt.figure(2)
        clf()
        axS = figS.add_subplot(1, 1, 1)

        axS.plot(prof[j].salt, prof[j].z, c="black")
        axS.plot(ctd[j].salt[gind], ctd[j].z[gind], c="red")
        axS.set_title(string("Station ", j, " (Salinity)"))
        axS.set_xlabel("Salinity")
        axS.set_ylabel("Depth (m)")
        axS.legend(["WOA", "NPP"])
        plt.savefig(string(figoutdir, "salt_sta_", string(j/10000)[end-1:end], "_", season[i], ".png"), dpi = 300)
    end
end

end #module
