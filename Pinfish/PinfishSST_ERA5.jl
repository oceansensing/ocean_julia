# This is a Julia translation of the Matlab code PinfishSST.m
# Donglai Gong 2019-06-07
#
using Statistics
using NCDatasets, Missings, Dates
using JLD2
using Plots

initflag = 1
loadflag = 1
analysisflag = 1
saveflag = 1
plotflag = 1

#yrrangeall = [yyyy for yyyy=2003:2023];
yrrangehipinfish = [2012; 2013; 2015];
yrrangelopinfish = [2014; 2016; 2017; 2018; 2019];
#yrrange1 = [2020];

year1 = 2016;
year2 = 2023;
syear1 = string(year1, pad=4);
syear2 = string(year2, pad=4);

month1 = 4;
month2 = 4;
smonth1 = string(month1, pad=2);
smonth2 = string(month2, pad=2);

#yearrange = yrrangeall; # change THIS to select different years!!!
#yearrange = year1:year2;
yearrange = year1:year2;
monthrange = month1:month2;

if initflag == 1 #| !exist('latMAB')
    datadir = "/Volumes/C2PO_Data/Data/ERA5/";
    #datapath = "https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2018/148/20180528090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc"
    datapath = datadir * "2022/ERA5_surface_hourly_202204.nc";
    ds = Dataset(datapath);
    lon = ds["longitude"];
    lat = ds["latitude"];
    t = ds["time"][:];
    lonarr = Float64.(collect(Missings.replace(lon[:],NaN)));
    latarr = Float64.(collect(Missings.replace(lat[:],NaN)));
    latind = findall( 44 .>= latarr .>= 33);
    lonind = findall(360.0-80.0 .<= lonarr .<= 360.0-69.0);
    latMAB = latarr[latind];
    lonMAB = lonarr[lonind];
end #if

#startind = [lonind[1],latind[1],1];
#countind = [size(lonind,1),size(latind,1),1];
#lonind2 = [ind for ind=startind[1]:(startind[1]+countind[1]-1)];
#latind2 = [ind for ind=startind[2]:(startind[2]+countind[2]-1)];
#tind = [ind for ind=startind[3]:(startind[3]+countind[3]-1)];

# load all the SST data and compute composite for MAB
if loadflag == 1
    sstMAB = repeat([0],size(latMAB,1),size(lonMAB,1));
    sstMAB1yr = sstMAB;
    sstMAByr = repeat([NaN],size(latMAB,1),size(lonMAB,1),size(yearrange,1));

    global navg = 0;

    #df = DateFormat("yyyymmddTHMS")

    for ii = 1:size(yearrange,1)
        yyyy = yearrange[ii];
        println(yyyy)
        navg1 = 0;

        for jj = 1:size(monthrange,1)
            global navg;
            global sstMAB;
            global sstMAB1yr;
            mm = monthrange[jj];
            println(mm)
            navg1 = navg1+1.0;
            navg = navg+1.0;

            datapath = datadir * string(yyyy) * "/" * "ERA5_surface_hourly_" * string(yyyy) * string(mm, pad=2) * ".nc";
            sstds = Dataset(datapath)["sst"];
            sstMAB1monthMiss = sstds[lonind[1]:lonind[end],latind[1]:latind[end],:];
            sstMAB1month = collect(Missings.replace(sstMAB1monthMiss,NaN));
            sstMAB1monthMean = mean(sstMAB1month, dims=3)[:,:,1];

            sstMAB1yr = sstMAB1yr*((navg1-1.0)/navg1) + sstMAB1monthMean/navg1;
            sstMAB = sstMAB*((navg-1.0)/navg) + sstMAB1monthMean/navg;
        end #for
        sstMAByr[:,:,ii] = sstMAB1yr;
    end #for
end #if

#latind2 = [ind for ind=startind[2]:(startind[2]+countind[2]-1)];
if analysisflag == 1
    LON = repeat(lonMAB',length(latMAB));
    LAT = repeat(reverse(latMAB,dims = 1),1,length(lonMAB));
    llind = findall((360.0-73.0 .>= LON .>= 360.0-75.0) .& (39 .>= LAT .>= 37));
    #nanmean(x) = mean(filter(!isnan,x)) # implementing nanmean similar to Matlab
    #nanmean(x,y) = mapslices(nanmean,x,dims = y) # implementing nanmean by dim of array
    sstMABts = mapslices(NaNMath.mean, sstMAByr[llind,:], dims=1)[:];
    #sstMABts = [];
    tyr = yearrange;
end #if

if saveflag == 1
    workdir = "/Users/gong/Research/Pinfish/";
    filepath = string(workdir, "sstMAB_", syear1, "-", syear2, "_", smonth1, "-", smonth2, ".jld2")
    jldsave(filepath; latMAB, lonMAB, sstMAB, sstMAByr, navg, tyr, sstMABts);
end #if

if plotflag == 1
    sstTSplot = Plots.plot(tyr, sstMABts  .- 272.15, framestyle=:box, title="ERA5 SST over southern MAB", label="Month: " * string(monthrange))
    Plots.savefig(sstTSplot, filepath[1:end-5] * ".html")
end #if
