% this script loads SST data and calculate composites from NASA's PO DAAC ghrsst server using MUR SST product.
% Example: https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2018/148/20180528090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc
% 4/16/2019: gong@vims.edu
% 10/11/2019: gong@vims.edu, updated to include 2019
% 10/14/2019: gong@vims.edu, updated to use nctoolbox

dbstop if error
close all
clear all

initflag = 1
loadflag = 1
csvflag = 1
analysisflag = 1
saveflag = 1

yrrangeall = [2003:2020];
yrrangehipinfish = [2012 2013 2015 2019];
yrrangelopinfish = [2014 2016 2017 2018];
yrrange1 = [2020];

%ydayStart = 1
%ydayEnd = 60
%ydayStart = 100
%ydayEnd = 160
ydayStart = 150
ydayEnd = 180

yday1 = ydayStart;
yday2 = ydayEnd;
syday1 = num2str(yday1+1000); syday1 = syday1(2:end);
syday2 = num2str(yday2+1000); syday2 = syday2(2:end);

ydayNlatest = 180

yearrange = yrrangeall % change THIS to select different years!!!
%yearrange = yrrangehipinfish
%yearrange = yrrangelopinfish
%yearrange = yrrange1

workdir = './';

if initflag == 1 | ~exist('latMAB')
    %% setting up datapath & netcdf paths
    datapath0 = 'https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2018/148/20180528090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc';
    ds = ncdataset(datapath0);

    %ncdisp(datapath);
    %datainfo = ncinfo(datapath);
    % determine the dimension of each variable
    %nt = datainfo.Variables(1).Dimensions.Length;
    %nlon = datainfo.Variables(3).Dimensions.Length;
    %nlat = datainfo.Variables(2).Dimensions.Length;
    %sst = repmat(NaN,[nlat,nlon,nt]);

    % set boundary for lon/lat for MAB
    %lon = double(ncread(datapath,'lon'));
    %lat = double(ncread(datapath,'lat'));
    lon = double(ds.data('lon'));
    lat = double(ds.data('lat'));
    latind = find(lat >= 33 & lat <= 41);
    lonind = find(lon <= -69 & lon >= -80);

    % extract and download data for MAB
    lonMAB = lon(lonind)';
    latMAB = lat(latind);

    %datapath = 'https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2018/148/*.nc';
    startind = [lonind(1),latind(1),1];
    countind = [length(lonind),length(latind),1];

    sstMAB = repmat(0,[length(latMAB),length(lonMAB)]);
    sstMAB1yr = sstMAB;
    sstMAByr = repmat(NaN,[length(latMAB),length(lonMAB),length(yearrange)]);

    t0 = datenum(2002,1,152);
    tN = datenum(2020,1,ydayNlatest);
    ydrange = [1:366];
    leapyr =@(yr)~rem(yr,400)|rem(yr,100)&~rem(yr,4);
    yday0 = 1;

    filenametail = ['090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'];

    yyyy1 = yearrange(1);
    yyyy2 = yearrange(end);

    for ii = 1:length(yearrange)
        yyyy = yearrange(ii)

        if leapyr(yyyy) == 1
            yday0 = 1;
            ydayN = 366;
        else
            yday0 = 1;
            ydayN = 365;
        end %if

        switch yyyy
        case 2002
            yday0 = 152;
            yday1 = yday0;
            yday2 = 160;
        case 2020
            yday1 = ydayStart;
            ydayN = ydayNlatest;
            if yday2 > ydayNlatest
                yday2 = ydayNlatest;
            end %if
        otherwise
            yday1 = ydayStart;
            yday2 = ydayEnd;
            %yday1 = yday0;
            %yday2 = ydayN;
        end %switch

        %yday1 = 90 % end of March
        %yday2 = 276 % beginning of October
        syday1 = num2str(yday1+1000); syday1 = syday1(2:end);
        syday2 = num2str(yday2+1000); syday2 = syday2(2:end);

        ydayarr = [yday1:yday2];

    end %for
end %if

% load all the SST data and compute composite for MAB
if loadflag == 1
    navg = 0;
    for ii = 1:length(yearrange)
        yyyy = yearrange(ii)
        navg1 = 0;
        sstMAB1yr = repmat(0,[length(latMAB),length(lonMAB)]);

        ydayarr = [yday1:yday2];

        sstMABi = repmat(NaN, [length(latMAB), length(lonMAB), length(ydayarr)]);

        for jj = 1:length(ydayarr)
            yday = ydayarr(jj)
            sstdata = [];
            filenamehead = datestr(datenum(yyyy,1,yday),30);
            filename = [filenamehead(1:8) filenametail]
            yday1000 = num2str(yday+1000);
            datapath = ['https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/' num2str(yyyy) '/' yday1000(2:end) '/' filename]
            navg = navg + 1
            navg1 = navg1 + 1;
            %sstdata0 = ncread(datapath,'analysed_sst',startind,countind)';

            startIdx = [1 startind(2) startind(1)];
            endIdx = [1 startind(2)+countind(2)-1 startind(1)+countind(1)-1];
            sstds = ncdataset(datapath);
            sstdata = squeeze(sstds.data('analysed_sst', startIdx, endIdx));
            sstds = [];
            sstMABi(:,:,jj) = sstdata;

            sstMAB1yr = sstMAB1yr*(1.0-1.0/navg1) + sstdata*(1.0/navg1);
            sstMAB = sstMAB*(1.0-1.0/navg) + sstdata*(1.0/navg);

        end %if
        sstMAByr(:,:,ii) = sstMAB1yr;
        save([workdir 'sstMABi_' num2str(yyyy) '_yd' syday1 '-' syday2 '.mat'],'sstMABi','lonMAB','latMAB','yyyy','ydayarr','navg1','-v7.3')
    end %for

    if saveflag == 1
        save([workdir 'sstMABm1_' num2str(yyyy1) '-' num2str(yyyy2) '_yd' syday1 '-' syday2 '.mat'],'sstMAB','sstMAByr','yearrange','lonMAB','latMAB','yday1','yday2','navg')
        %save([workdir 'sstMAB_' num2str(yyyy1) '-' num2str(yyyy2) '.mat'],'sstMAB','lonMAB','latMAB','navg')
    end %if

end %if loadflag

if csvflag == 1
    ydind = [yday1:30:yday2]-yday1+1;
    ydd1 = ydind(1:end-1); ydd1(1) = ydind(1);
    ydd2 = ydind(2:end)-1; ydd2(end) = ydind(end);
    ydd1start = ydd1+yday1-1;

    csvmatrixSlopeSea = repmat(NaN,[length(yrrangeall), length(ydd1)]);
    csvmatrixChesapeake = repmat(NaN,[length(yrrangeall), length(ydd1)]);

    % calculating the monthly average SST for each year for two regions in the MAB
    for ii = 1:length(yrrangeall)
        ii
        yyyy = yearrange(ii)
        load([workdir 'sstMABi_' num2str(yyyy) '_yd' syday1 '-' syday2 '.mat']);

        [LON,LAT] = meshgrid(lonMAB,latMAB);
        sstMABtsSlopeSea = [];
        sstMABtsChesapeake = [];
        llindSlopeSea = find(LON >= -74.5 & LON <= -73 & LAT >= 36.5 & LAT <= 39);
        llindChesapeake = find(LON >= -77 & LON <= -75 & LAT >= 36.5 & LAT <= 39);
        for jj = 1:length(ydd1)
            sstMAB1month = sstMABi(:,:,[ydd1(jj):ydd2(jj)])-273;
            csvmatrixSlopeSea(ii,jj) = nanmean(sstMAB1month(llindSlopeSea));
            csvmatrixChesapeake(ii,jj) = nanmean(sstMAB1month(llindChesapeake));
        end %for
    end %for

    %formatting matrix for writing to CSV files
    csvmatrixSlopeSea = [yrrangeall' , csvmatrixSlopeSea];
    csvmatrixSlopeSea = [[NaN ydd1start] ; csvmatrixSlopeSea];

    csvmatrixChesapeake = [yrrangeall' , csvmatrixChesapeake];
    csvmatrixChesapeake = [[NaN ydd1start] ; csvmatrixChesapeake];

    csvwrite(['sstMABcsv_SlopeSea_yd' syday1 '-' syday2 '.csv'],csvmatrixSlopeSea);
    csvwrite(['sstMABcsv_Chesapeake_yd' syday1 '-' syday2 '.csv'],csvmatrixChesapeake);

end %if

if analysisflag == 1
    load(['sstMABm1_2003-2020_yd' syday1 '-' syday2 '.mat']);
    syday1 = num2str(yday1+1000); syday1 = syday1(2:end);
    syday2 = num2str(yday2+1000); syday2 = syday2(2:end);

    %sstMAB = repmat(0,[length(latMAB),length(lonMAB)]);
    nnanind = find(~isnan(sstMAB));
    sstMAB(nnanind) = 0.0;
    sstMAB1yr = sstMAB;
    sstMAByr = repmat(NaN,[length(latMAB),length(lonMAB),length(yearrange)]);

    navg = 0;
    for ii = 1:length(yearrange)
        ii
        yyyy = yearrange(ii)
        load([workdir 'sstMABi_' num2str(yyyy) '_yd' syday1 '-' syday2 '.mat']);
        navg = navg + 1;
        sstMAB1yr = nanmean(sstMABi,3);
        pcolor(sstMAB1yr); shading flat; colormap(jet(100));
        sstMAByr(:,:,ii) = sstMAB1yr;
        sstMAB = sstMAB*(1.0-1.0/navg) + sstMAB1yr*(1.0/navg);
        %pause;
    end %for
 %    pcolor(lonMAB,latMAB,sstMAByr(:,:,1)); shading flat; colormap(jet(100)); colorbar;

    [LON,LAT] = meshgrid(lonMAB,latMAB);
    sstMABts = [];
    sstMABii = [];
    llind = find(LON >= -74.5 & LON <= -73 & LAT >= 37.5 & LAT <= 39);
    for ii = 1:size(sstMAByr,3)
        sstMABii = sstMAByr(:,:,ii);
        sstMABts(ii) = nanmean(sstMABii(llind));
    end %for
    tyr = yearrange;

    if saveflag == 1
        save([workdir 'sstMABm2_' num2str(yyyy1) '-' num2str(yyyy2) '_yd' syday1 '-' syday2 '.mat'],'sstMAB','sstMAByr','sstMABts','tyr','yearrange','lonMAB','latMAB','yday1','yday2','navg','-v7.3')
        %save([workdir 'sstMAB_' num2str(yyyy1) '-' num2str(yyyy2) '.mat'],'sstMAB','lonMAB','latMAB','navg')
    end %if
end %if
