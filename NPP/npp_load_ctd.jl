module npp_load_ctd
# 2019-10-01 Donglai Gong
# 2019-11-08 DG, ported over more MATLAB version's functions
# 2019-11-13 DG, cleaned up the code some
# 2019-11-19 DG, further cleaned up the code and make save/load in jld2 format work.
#
# Note: you should not need to call this module directly, use "run_npp_load_ctd.jl" instead

# load the needed packages
using CSV, Glob, FileIO, JLD2 #, DataFrames, Serialization, FileIO, JLD
using Missings, Dates, GSW, Statistics, DataFrames
#using Random, Distributions
using PyCall, WebIO
#np = pyimport("numpy");
sklen = pyimport("sklearn.ensemble");

# add the work directory into LOAD_PATH if it's not already done
if (pwd() in LOAD_PATH) == false
    push!(LOAD_PATH, pwd());
end

# importing and exporting all the relevant variables needed for this analysis
import C2PO: gc_distance, rad2deg, deg2rad, histc, meshgrid, nan, findNaNmin, findNaNmax # this module stores additional function needed for analysis (e.g. gc_distance)
import NPP_types: CTD, Expedition # this module stores all the needed data Types (i.e. data structure)
export CTD, ctd, Expedition, npp, datapath, load_ctd_cnv, load_npp, reloadflag, saveflag # lists the variables that can be accessed outside the module

# this file defines the flags controlling loading and saving
include("npp_flags.jl")

# set the reloadflag and saveflag if they are not defined already
if (@isdefined reloadflag) == false
    global reloadflag = 1
end

if (@isdefined saveflag) == false
    global saveflag = 1
end

println(string("reloadflag: ", reloadflag))
println(string("saveflag: ", saveflag))

# load the file that defines the data and code directory paths to be used
include("npp_path_setup.jl")

# calculate the number of casts there are and store in the ncasts variable.
const ncasts = size(ctddatapath,1);
global const puplim = 5.0; # shallowest depth to analyze data;
global const pdnlim = 1000.0; # deepest depth to anayze data

# function for loading the CTD CNV data files
function load_ctd_cnv(ctddatapath, CTD)
    ctdheader = ["timeJ", "timeY", "latitude", "longitude", "prDM", "depSM", "t090C", "t190C", "c0S_m", "c1S_m", "sal00", "sal11", "svCM", "svCM1", "density00", "sigma00", "density11", "sigma11", "sbeox0V", "sbeox0Mg_L", "sbeox0ML_L", "sbeox0PS", "turbWETntu0", "flECO_AFL", "altM", "sfdSM", "accM", "dz_dtM", "bpos", "pumps", "flag"];
    yday = Array{Float64};
    lon = Array{Float64};
    lat = Array{Float64};
    p = Array{Float64};
    z = Array{Float64};
    temp0 = Array{Float64};
    temp1 = Array{Float64};
    temp = Array{Float64};
    ptmp = Array{Float64};
    ctmp = Array{Float64};
    salt0 = Array{Float64};
    salt1 = Array{Float64};
    salt = Array{Float64};
    saltA = Array{Float64};
    rho = Array{Float64};
    sigma0 = Array{Float64};
    oxy = Array{Float64};
    oxysat = Array{Float64};
    fluor = Array{Float64};
    casttype = Array{Int64};
    dataflag = Array{Int64};

    ctd = CTD[];
    maxrow = 1;

    for ii = 1:ncasts
        println(ctddatapath[ii])
        ctdcnv = CSV.read(ctddatapath[ii], header=ctdheader, skipto=1330, delim=' ', ignoreemptylines=true, ignorerepeated=true, comment="#", missingstring="NaN");
        maxrow = max(size(ctdcnv,1),maxrow)
        yday = Array(ctdcnv.timeJ);

        lon = Array(ctdcnv.longitude);
        lat = Array(ctdcnv.latitude);

        p = Array(ctdcnv.prDM);
        zctd = -Array(ctdcnv.depSM);
        #z = map(GSW.gsw_z_from_p, p.+10.1325, lat);
        z = GSW.gsw_z_from_p.(p.+10.1325, lat); # strange bug with the first element of z, need to manually reset it.
        z[1] = z[2];

        bind = findall((p .> pdnlim) .| (p .< -10));
        z[bind] .= NaN;
        p[bind] .= NaN;

        salt0 = Array(ctdcnv.sal00);
        salt1 = Array(ctdcnv.sal11);
        salt = (salt0 .+ salt1) ./2;
        saltA = GSW.gsw_sa_from_sp.(salt, p, lon, lat);

        temp0 = Array(ctdcnv.t090C);
        temp1 = Array(ctdcnv.t190C);
        temp = (temp0 .+ temp1) ./2;
        ctmp = GSW.gsw_ct_from_t.(saltA, temp, p);
        ptmp = GSW.gsw_pt_from_t.(saltA, temp, p, zeros(size(p,1)));
        rho = GSW.gsw_rho.(saltA, ctmp, p);
        sigma0 = GSW.gsw_sigma0.(saltA, ctmp);
        oxy = Array(ctdcnv.sbeox0Mg_L);
        oxysat = Array(ctdcnv.sbeox0PS);
        fluor = Array(ctdcnv.flECO_AFL);

        (zmin,zmind) = findNaNmax(p);
        casttype = Array{Int64,1}(undef,length(p));
        casttype[1:zmind] .= -1; # downcast
        casttype[zmind+1:end] .= +1; # upcast
        #pind = findall((casttype .< 0) .& (z .< -3));
        #pind = findall((casttype .< 0) .& (p .>= 5));
        pind = findall(p .>= puplim);
        dataflag = [-1 for ii in 1:length(p)];
        #dataflag[pind] .= 1;

        push!(ctd, CTD(yday,lon,lat,p, z, temp, ptmp, ctmp, salt, saltA, rho, sigma0, oxy, oxysat, fluor,casttype,pind,dataflag));
    end
    clean_ctd!(ctd);
    return ctd
end

# function for loading the CTD CSV data files (no need to run this and cnv)
function load_ctd_csv(ctddatapathcsv, CTD)
    maxrow = 1;
    ncasts = size(ctddatapathcsv,1);
    yday = Array{Float64};
    lon = Array{Float64};
    lat = Array{Float64};
    p = Array{Float64};
    z = Array{Float64};
    temp = Array{Float64};
    ptmp = Array{Float64};
    ctmp = Array{Float64};
    salt = Array{Float64};
    saltA = Array{Float64};
    rho = Array{Float64};
    sigma0 = Array{Float64};
    oxy = Array{Float64};
    oxysat = Array{Float64};
    fluor = Array{Float64};
    dataflag = Array{Float64};

    # initialize ctd data array
    ctd = CTD[]
    for ii = 1:ncasts
        maxrow
        ctdcsv = CSV.read(ctddatapathcsv[ii], skipto=2, comment="#", missingstring="NaN")
        maxrow = max(size(ctdcsv,1),maxrow)
        yday = Float64.(collect(Missings.replace(ctdcsv[:,1],NaN)));
        lon = Float64.(collect(Missings.replace(ctdcsv[:,2],NaN)));
        lat = Float64.(collect(Missings.replace(ctdcsv[:,3],NaN)));
        p = Float64.(collect(Missings.replace(ctdcsv[:,4],NaN)));
        z = Float64.(collect(Missings.replace(ctdcsv[:,5],NaN)));
        temp = Float64.(collect(Missings.replace(ctdcsv[:,6],NaN)));
        ptmp = Float64.(collect(Missings.replace(ctdcsv[:,9],NaN)));
        ctmp = Float64.(collect(Missings.replace(ctdcsv[:,10],NaN)));
        salt = Float64.(collect(Missings.replace(ctdcsv[:,8],NaN)));
        saltA = Float64.(collect(Missings.replace(ctdcsv[:,12],NaN)));
        rho = Float64.(collect(Missings.replace(ctdcsv[:,13],NaN)));
        sigma0 = Float64.(collect(Missings.replace(ctdcsv[:,14],NaN)));
        oxy = Float64.(collect(Missings.replace(ctdcsv[:,17],NaN)));
        oxysat = Float64.(collect(Missings.replace(ctdcsv[:,19],NaN)));
        fluor = Float64.(collect(Missings.replace(ctdcsv[:,20],NaN)));

        bind = findall((p .> pdnlim) .| (p .< -10));
        z[bind] .= NaN;
        p[bind] .= NaN;

        (zmin,zmind) = findNaNmin(z);
        casttype = Array{Int64,1}(undef,length(z));
        casttype[1:zmind] .= -1; # downcast
        casttype[zmind+1:end] .= +1; # upcast
        #pind = findall((casttype .< 0) .& (z .< -3));
        pind = findall(p .>= puplim);
        dataflag = [-1 for ii in 1:length(z)];
        #dataflag[pind] .= 1;
        #println(string(ii) * " " * string(maxrow))

        # adding data to the ctd data structure of composite type CTD
        push!(ctd, CTD(yday,lon,lat,p, z, temp, ptmp, ctmp, salt, saltA, rho, sigma0, oxy, oxysat, fluor, casttype, pind, dataflag));
    end
    clean_ctd!(ctd);
    return ctd
end

# load_npp loads the cruise meta data. output as npp data type as defined above
function load_npp(ctd)
    ncasts = length(ctd);
    project = "NPP"
    chiefsci = "Brice Loose"
    sections = [[1 2 3 4], [5], [6 7 8 9 10 11 12 13 14 15 16], [17 18 19 20 21 22], [23], [24 25 26 27 28 29], [31 32 30 33 34 35 36 37 38 47 48 49 50], [39 40 41 42 43 44 45 46], [51 52]];
    sectnames = ["JonesSound", "PondInlet", "LancasterSound", "WellingtonChannel", "MelvilleSound", "PeelSound", "BarrowStraitEast", "PrinceRegentInlet", "CrokerBay"];

    # calculate the average lon/lat locations for each CTD station
    mlon = Array{Float64,1}(undef,ncasts);
    mlat = Array{Float64,1}(undef,ncasts);
    for ii = 1:ncasts
        mlon[ii] = mean(ctd[ii].lon);
        mlat[ii] = mean(ctd[ii].lat);
    end #for

    # calculate the along transect distance for each transect
    stransect = [];
    for ii = 1:length(sections)
        sind = sections[ii];
        if length(sind) > 1
            push!(stransect,[0;cumsum(gc_distance.(mlat[sind[1:end-1]],mlon[sind[1:end-1]],mlat[sind[2:end]],mlon[sind[2:end]]))]);
        else
            push!(stransect,[0]);
        end
    end #for

    # save project and transect info in npp data structure of type Expedition
    npp = Expedition(project,chiefsci,sections,sectnames,Date(2019,7,15),Date(2019,8,5),stransect); #This is the constructor for the structure Expedition. npp is immuatable.
    return npp
end

# calculate dataflag for CTD profiles using Isolation Forest algorithm for problematic casts, otherwise whitelist non problematic downcast. dataflag of +1 is good, -1 is bad.
# clean up CTD temp and salinity data by flagging those outside realistic range
function clean_ctd!(ctd)
    problemcasts = [5 10 15 16 20 22 23 25 27 28 29 30 31 32 34 35 36 38 41 42 43 45 47 49 51]
    noprobcasts = setdiff(collect(1:length(ctd)), vec(problemcasts))

    if !isempty(problemcasts)
        for ii in problemcasts
            #ii = 28
            println(string("CTD cast #: ", ii))
            pindraw = ctd[ii].pind;
            #gind = findall((ctd[ii].ctmp[pindraw] .< 3) .& (ctd[ii].ctmp[pindraw] .> -1.95) .& (ctd[ii].saltA[pindraw] .< 35) .& (ctd[ii].saltA[pindraw] .> 26) .& (ctd[ii].z[pindraw] .> -800));

            gindShallowExt = findall((140 .>= ctd[ii].p[pindraw] .>= puplim-1.0) .& (ctd[ii].casttype[pindraw] .< 0));
            gindDeepExt = findall((pdnlim .> ctd[ii].p[pindraw] .> 70) .& (ctd[ii].casttype[pindraw] .< 0));

            gindShallow = findall((90 .>= ctd[ii].p[pindraw] .>= puplim) .& (ctd[ii].casttype[pindraw] .< 0));
            gindDeep = findall((pdnlim .> ctd[ii].p[pindraw] .> 90) .& (ctd[ii].casttype[pindraw] .< 0));

            #gindShallowExt = findall(140 .>= ctd[ii].p[pindraw] .>= 3);
            #gindDeepExt = findall(ctd[ii].p[pindraw] .> 70);

            #gindShallow = findall(90 .>= ctd[ii].p[pindraw] .>= 5 );
            #gindDeep = findall(ctd[ii].p[pindraw] .> 90);

            pindShallowExt = pindraw[gindShallowExt];
            pindDeepExt = pindraw[gindDeepExt];

            pindShallow = pindraw[gindShallow];
            pindDeep = pindraw[gindDeep];

            pind = [pindShallow; pindDeep];

            ctdShallowExt_data = [ctd[ii].p[pindShallowExt] ctd[ii].saltA[pindShallowExt]];
            ctdDeepExt_data = [ctd[ii].p[pindDeepExt] ctd[ii].saltA[pindDeepExt]];

            ctdShallow_data = [ctd[ii].p[pindShallow] ctd[ii].saltA[pindShallow]];
            ctdDeep_data = [ctd[ii].p[pindDeep] ctd[ii].saltA[pindDeep]];

            #ctd_train = [ctd[ii].p[pind[1:4:end]] ctd[ii].ctmp[pind[1:4:end]] ctd[ii].saltA[pind[1:4:end]]];
            #ctd_data = [ctd[ii].p[pind] ctd[ii].ctmp[pind] ctd[ii].saltA[pind]];

            clfShallow = sklen.IsolationForest(behaviour="new", max_samples=100000, contamination=0.0008)
            clfDeep = sklen.IsolationForest(behaviour="new", max_samples=100000, contamination=0.002, bootstrap=true)

            clfShallow.fit(ctdShallowExt_data);
            clfDeep.fit(ctdDeepExt_data);

            ctdShallow_pred_data = clfShallow.predict(ctdShallow_data);
            ctdDeep_pred_data = clfDeep.predict(ctdDeep_data);

            ctd_pred_data = [ctdShallow_pred_data; ctdDeep_pred_data];
            println(string("# of data total: ", length(pindraw)))
            println(string("# of data evals: ", length(ctd_pred_data)))
            #ctd[ii].dataflag[pindraw] .= ctd_pred_data; # 1 is good data, -1 is bad data
            ctd[ii].dataflag[pind] .= ctd_pred_data; # 1 is good data, -1 is bad data
        end
    end

    if !isempty(noprobcasts)
        for ii in noprobcasts
            pindraw = ctd[ii].pind;
            gind = findall((pdnlim .>= ctd[ii].p[pindraw] .>= puplim) .& (ctd[ii].casttype[pindraw] .< 0));
            ctd[ii].dataflag[pindraw[gind]] .= 1;
        end
    end

    # remove "bad" data points from all the CTD casts
    for ii in 1:length(ctd)
        bind = findall( (ctd[ii].p .< 0.0) .| (ctd[ii].salt .> 38.0) .| (ctd[ii].salt .< 20.0) .| (ctd[ii].temp .> 16.0) .| (ctd[ii].temp .< -2) .| ((ctd[ii].salt .< 31.0) .& (ctd[ii].p .> 50)) );
        ctd[ii].dataflag[bind] .= -1;
    end
end #clean_ctd!

# load ctd data and npp metadata, return both
if reloadflag == 1
    #ctd = load_ctd_csv(ctddatacsvpath,CTD); # uncomment this line if you want to load the CSV version
    ctd = load_ctd_cnv(ctddatapath,CTD);
    npp = load_npp(ctd);

    lon = Array{Float64,1}(undef,length(ctd));
    lat = Array{Float64,1}(undef,length(ctd));
    yday = Array{Float64,1}(undef,length(ctd));

    for i = 1:length(ctd)
        lon[i] = mean(ctd[i].lon);
        lat[i] = mean(ctd[i].lat);
        yday[i] = mean(ctd[i].yday);
    end
    df = DataFrame(yday = yday, lon = lon, lat = lat);
    CSV.write(ctdlocationpath, df)
    println("npp and ctd are reloaded from CNV data files.")
else
    project = load(filepath)
    ctd = project["ctd"];
    npp = project["npp"];
    println("npp and ctd are loaded in Julia.")
end

# save ctd and npp data and data types
if saveflag == 1 & reloadflag == 1
    #@time save(filepath, "ctd", ctd, "npp", npp)
    @time @save filepath ctd npp
    println("npp and ctd are saved to disk.")
end #if

end
