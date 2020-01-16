# This is a Julia translation of the Matlab code PinfishSST.m
# Donglai Gong 2019-06-07
#
using Statistics
using NCDatasets, Missings, Dates
using FileIO
using PlotlyJS

initflag = 1
loadflag = 1
analysisflag = 1
saveflag = 1

yrrangeall = [yyyy for yyyy=2003:2018];
yrrangehipinfish = [2012 2013 2015];
yrrangelopinfish = [2014 2016 2017 2018];
yrrange1 = [2016];

ydayNlatest = 157
yearrange = yrrangeall; # change THIS to select different years!!!

if initflag == 1 #| !exist('latMAB')
    datapath = "https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2018/148/20180528090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc"
    ds = Dataset(datapath);
    lon = ds["lon"];
    lat = ds["lat"];
    lonarr = Float64.(collect(Missings.replace(lon[:],NaN)));
    latarr = Float64.(collect(Missings.replace(lat[:],NaN)));
    latind = findall( 44 .>= latarr .>= 33);
    lonind = findall(-80 .<= lonarr .<= -69);
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

    t0 = Date(2002,1,1)+Day(152-1)
    tN = Date(2019,1,1)+Day(157-1) #dayofyear(Date(2019,5,5))

    yday0 = 1;
    global navg = 0;

    filenametail = "090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc";

    yyyy1 = yearrange[1];
    yyyy2 = yearrange[end];

    #df = DateFormat("yyyymmddTHMS")

    for ii = 1:size(yearrange,1)
        yyyy = yearrange[ii];
        println(yyyy)
        navg1 = 0;
        sstMAB1yr = repeat([0],size(latMAB,1),size(lonMAB,1));

        if Dates.isleapyear(Date(yyyy))
            yday0 = 1;
            ydayN = 366;
        else
            yday0 = 1;
            ydayN = 365;
        end #if

        yday1 = 120;
        yday2 = 160;

        if yyyy == 2002
            yday0 = 152;
            yday1 = yday0;
        elseif yyyy == 2019
            ydayN = ydayNlatest;
            yday2 = ydayN;
        else

        end

        yday1 = 90
        yday2 = ydayNlatest

        syday1000 = string(yday1+1000);
        syday2000 = string(yday2+1000);
        ydays = [yday for yday in yday1:yday2];
        for jj = 1:length(ydays)
            global navg
            global sstMAB
            yday = ydays[jj];
            println(yday)
            t = DateTime(yyyy,1,1)+Day(yday-1);
            st = string(t);
            filenamehead = string(st[1:4], st[6:7], st[9:10]);
            filename = string(filenamehead, filenametail);
            yday000 = string(yday+1000);
            datapath = string("https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/", yyyy, "/", yday000[2:end], "/", filename);
            
            navg = navg + 1;
            navg1 = navg1 + 1;

            sstds = Dataset(datapath)["analysed_sst"];
            sstMABi = sstds[lonind[1]:lonind[end],latind[1]:latind[end],1];
            sstMABimx = collect(Missings.replace(sstMABi,NaN));

            sstMAB1yr = sstMAB1yr*(1.0-1.0/navg1) + sstMABimx*(1.0/navg1);
            sstMAB = sstMAB*(1.0-1.0/navg) + sstMABimx*(1.0/navg);
        end #for
        sstMAByr[:,:,ii] = sstMAB1yr;
    end #for
end #if

#latind2 = [ind for ind=startind[2]:(startind[2]+countind[2]-1)];
if analysisflag == 1
    LON = repeat(lonMAB',length(latMAB));
    LAT = repeat(reverse(latMAB,dims = 1),1,length(lonMAB));
    llind = findall((-73 .>= LON .>= -74.5) .& (39 .>= LAT .>= 37.5));
    nanmean(x) = mean(filter(!isnan,x)) # implementing nanmean similar to Matlab
    nanmean(x,y) = mapslices(nanmean,x,dims = y) # implementing nanmean by dim of array
    #sstMABts = nanmean(sstMAByr[llind],)
    sstMABts = [];
    tyr = yearrange;
end #if

if saveflag == 1
    workdir = "/Users/gong/Documents/Research/Projects/Pinfish/Pinfish/";
    filepath = string(workdir, "sstMAB_", yyyy1, "-", yyyy2, "_yd", syday1000[2:end], "-", syday2000[2:end], ".jld2")
    @save filepath latMAB lonMAB sstMAB sstMAByr navg tyr sstMABts
end #if    

if plotflag == 1
    function sstmap()
        sst = heatmap(
            x=lonMAB,
            y=latMAB,
            z=sstMAB1yr
        )
        plot(sst)
    end
    sstmap()
end #if