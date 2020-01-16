module WOA

using NCDatasets, Glob, CSV, Missings

if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

#import C2PO: gc_distance, rad2deg, deg2rad, histc, meshgrid, nan, findNaNmin, findNaNmax # this module stores additional function needed for analysis (e.g. gc_distance)

export npp_load_woa, load_woa, season, Profile, WOAdata

include("npp_path_setup.jl")

mutable struct Profile
    z::Array{Float64}
    temp::Array{Float64}
    salt::Array{Float64}
end

mutable struct WOAdata
    lat::Array{Float64}
    lon::Array{Float64}
    z::Array{Float64}
    temp::Array{Float64}
    salt::Array{Float64}
end

# 1: winter, 2: spring, 3: summer, 4: fall
const season = Dict([(1,"Winter"), (2,"Spring"), (3,"Summer"), (4,"Fall")]);
const seasontempfile = Dict([(1,"woa18_decav81B0_t13_04.nc"), (2,"woa18_decav81B0_t14_04.nc"), (3,"woa18_decav81B0_t15_04.nc"), (4,"woa18_decav81B0_t16_04.nc")]);
const seasonsaltfile = Dict([(1,"woa18_decav81B0_s13_04.nc"), (2,"woa18_decav81B0_s14_04.nc"), (3,"woa18_decav81B0_s15_04.nc"), (4,"woa18_decav81B0_s16_04.nc")]);

# this function loads the global WOA dataset, not necessary (overkill) for comparison with NPP CTD data.
function load_woa(i::Int64 = 3, woadir::String=string(ENV["HOME"], "/Research/data/WOA18/"))
    #seasonlocal = season;
    println(season[i])

    temproot = "https://data.nodc.noaa.gov/thredds/dodsC/ncei/woa/temperature/decav81B0/0.25/";
    saltroot = "https://data.nodc.noaa.gov/thredds/dodsC/ncei/woa/salinity/decav81B0/0.25/";

    #woatemppath = Glob.glob("woa*81B0_t*.nc", woadir);
    #woasaltpath = Glob.glob("woa*81B0_s*.nc", woadir);
    woatemppath = temproot * seasontempfile[i];
    woasaltpath = saltroot * seasonsaltfile[i];

    # load WOA18 location data for the specified season
    dstemp = Dataset(woatemppath);
    dssalt = Dataset(woasaltpath);
    latnc = dstemp["lat"];
    lonnc = dstemp["lon"];
    depthnc = dstemp["depth"];
    lon = Float64.(collect(Missings.replace(lonnc[:],NaN)));
    lat = Float64.(collect(Missings.replace(latnc[:],NaN)));
    z = -Float64.(collect(Missings.replace(depthnc[:],NaN)));
    #(LON,LAT) = meshgrid(lon,lat);

    tempnc = dstemp["t_an"];
    saltnc = dssalt["s_an"];
    temp = Float64.(collect(Missings.replace(tempnc[:],NaN)));
    salt = Float64.(collect(Missings.replace(saltnc[:],NaN)));
    woa = WOAdata(lat, lon, z, temp, salt);
    #return lat, lon, z, temp, salt;
    return woa;
end

# this function load the WOA data only at the closest WOA location for each NPP CTD station (requires internet connection)
function npp_load_woa(i::Int64 = 3, ctdlocationpath::String = ctdlocationpath)
    # load WOA data nearest to NPP CTD station for a specific season
    # i = 1: winter, 2: spring, 3: summer, 4: fall

    #woadir = string(ENV["HOME"], "/Research/data/WOA18/");
    #ctdlocationpath = string(ENV["HOME"], "/Research/NPP/data/CTD/NPP_ctd_locations.csv")

    #woa = load_woa(i,woadir); # only use this if want to load the global WOA dataset.
    temproot = "https://data.nodc.noaa.gov/thredds/dodsC/ncei/woa/temperature/decav81B0/0.25/";
    saltroot = "https://data.nodc.noaa.gov/thredds/dodsC/ncei/woa/salinity/decav81B0/0.25/";

    #woatemppath = Glob.glob("woa*81B0_t*.nc", woadir); # use this if a version of the WOA data is downloaded locally
    #woasaltpath = Glob.glob("woa*81B0_s*.nc", woadir);
    woatemppath = temproot * seasontempfile[i];
    woasaltpath = saltroot * seasonsaltfile[i];

    # load WOA18 location data for the specified season
    dstemp = Dataset(woatemppath);
    dssalt = Dataset(woasaltpath);
    latnc = dstemp["lat"];
    lonnc = dstemp["lon"];
    depthnc = dstemp["depth"];
    lon = Float64.(collect(Missings.replace(lonnc[:],NaN)));
    lat = Float64.(collect(Missings.replace(latnc[:],NaN)));
    z = -Float64.(collect(Missings.replace(depthnc[:],NaN)));
    #(LON,LAT) = meshgrid(lon,lat);

    tempnc = dstemp["t_an"];
    saltnc = dssalt["s_an"];
    #temp = Float64.(collect(Missings.replace(tempnc[:],NaN)));
    #salt = Float64.(collect(Missings.replace(saltnc[:],NaN)));
    #woa = WOAdata(lat, lon, z, temp, salt);

    #woa = load_woa(i, woadir);

    # location NPP CTD location data
    ctdloc = CSV.read(ctdlocationpath); # reading in the NPP CTD location data
    lonctd = convert(Vector,ctdloc.lon); # converting from dataframe to vector
    latctd = convert(Vector,ctdloc.lat);

    prof = Profile[];

    for j = 1:length(lonctd)
        println("Loading Station: " * string(j))
        #latind = findmin(abs.(woa.lat .- latctd[j]))[2];
        #lonind = findmin(abs.(woa.lon .- lonctd[j]))[2];;
        #push!(prof,Profile(woa.z, woa.temp[lonind,latind,:,1], woa.salt[lonind,latind,:,1]));
        latind = findmin(abs.(lat .- latctd[j]))[2];
        lonind = findmin(abs.(lon .- lonctd[j]))[2];

        temp = collect(Missings.replace(tempnc[lonind,latind,:,1],NaN));
        salt = collect(Missings.replace(saltnc[lonind,latind,:,1],NaN));
        push!(prof, Profile(z, temp, salt));
    end
    return prof
end

end #module
