% use this script to plot the SST results from running PinfishSST.m
% please run PinfishSST.m first!
%
%2019-04-18: gong@vims.edu
%
analysisflag = 1

yrrangeall = [2003:2020];
yrrangehipinfish = [2012 2013 2015 2019];
yrrangelopinfish = [2014 2016 2017 2018];
[tmp,yrhiind] = intersect(yrrangeall,yrrangehipinfish);
[tmp,yrloind] = intersect(yrrangeall,yrrangelopinfish);
yyyy = 2020;

%sstlop5 = load('sstMAB_2014-2018_yd100-160.mat');
%ssthip5 = load('sstMAB_2012-2019_yd100-160.mat');
%sstnow2015 = load('sstMAB_2015-2015_yd100-160.mat');
%sstnow2016 = load('sstMAB_2016-2016_yd100-160.mat');
%sstnow2019 = load('sstMAB_2019-2019_yd100-160.mat');
load('sstMABm2_2003-2019_yd100-160.mat');

sstMABhi = nanmean(sstMAByr(:,:,yrhiind),3);
sstMABlo = nanmean(sstMAByr(:,:,yrloind),3);

figfmt = 'png' %png or epsc

% 1 = avg, 2 = anom, 4 = hilo, 8 = yearlyanom, 16 = now, 32 = year1anom
plotflag = 4

if analysisflag == 1
    [LON,LAT] = meshgrid(lonMAB,latMAB);
    sstMABi = [];
    sstMABtsSlopeSea = [];
    sstMABtsChesapeake = [];
    llindSlopeSea = find(LON >= -74.5 & LON <= -73 & LAT >= 36.5 & LAT <= 39);
    llindChesapeake = find(LON >= -77 & LON <= -75 & LAT >= 36.5 & LAT <= 39);
    for ii = 1:size(sstMAByr,3)
        sstMABi = sstMAByr(:,:,ii);
        sstMABtsSlopeSea(ii) = nanmean(sstMABi(llindSlopeSea));
        sstMABtsChesapeake(ii) = nanmean(sstMABi(llindChesapeake));
    end %for
    sstMABtsSlopeSea = sstMABtsSlopeSea - 273;
    sstMABtsChesapeake = sstMABtsChesapeake - 273;
    tyr = yrrangeall;
end %if

% time series plot
if ismember(plotflag,[0])
    fTS = figure('unit','inches');
    %hp = pcolor(lonMAB,latMAB,sstMAB-273); shading flat;
    hp = plot(tyr,sstMABtsSlopeSea,'*-')
    set(hp,'linewidth',1)
    set(gca,'box','on','tickdir','out','xgrid','on','ygrid','on')
    ht = title('April & May MAB Slope Sea SST (2003-2019)');
    hx = xlabel('Time');
    hy = ylabel('Sea Surface Temperature');
    set(gca,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    print('figPinfishSST_timeseriesSlopeSea',['-d' figfmt],'-r300');
    close(fTS)

    fTS = figure('unit','inches');
    %hp = pcolor(lonMAB,latMAB,sstMAB-273); shading flat;
    hp = plot(tyr,sstMABtsChesapeake,'*-')
    set(hp,'linewidth',1)
    set(gca,'box','on','tickdir','out','xgrid','on','ygrid','on')
    ht = title('April & May Chesapeake Bay SST (2003-2019)');
    hx = xlabel('Time');
    hy = ylabel('Sea Surface Temperature');
    set(gca,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    print('figPinfishSST_timeseriesChesapeake',['-d' figfmt],'-r300');
    close(fTS)

end %if

% avg
if ismember(plotflag,[1 3 6 9])
    favg5 = figure('unit','inches');
    hp = pcolor(lonMAB,latMAB,sstMAB-273); shading flat;
    set(gca,'box','on','tickdir','out')
    caxis([8 27]); colorbar;
    colormap(jet((round((27-8)/0.1))));
    ht = title('April & May MAB Slope Sea SST (2003-2019)');
    hx = xlabel('Longitude');
    hy = ylabel('Latitude');
    set(gca,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    print('figPinfishSST_avg5',['-d' figfmt],'-r300');
    close(favg5)
end %if

% now
if ismember(plotflag,[16 17 18 20 24])
    [tmp,yr1ind] = intersect(yrrangeall,yyyy);
    fnow = figure('unit','inches')
    pcolor(lonMAB,latMAB,sstMAByr(:,:,yr1ind) - sstMAB); shading flat;
    set(gca,'box','on','tickdir','out')
    caxis([-6 6]); colorbar;
    colormap(jet(10/0.1));
    ht = title(['Spring ' num2str(yyyy) ' SST Anomaly (yday 100-160)']);
    hx = xlabel('Longitude');
    hy = ylabel('Latitude');
    set(gca,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    print(['figPinfishSST_now' num2str(yyyy)],['-d' figfmt],'-r300');
    close(fnow)
end %if

% anom
if ismember(plotflag,[2 3 6 10])
    fhianom = figure('unit','inches');
    pcolor(lonMAB,latMAB,sstMABhi-sstMAB); shading flat;
    set(gca,'box','on','tickdir','out')
    caxis([-3 3]); colorbar;
    colormap(jet(6/0.1));
    ht = title('SST anomaly of high Pinfish years');
    hx = xlabel('Longitude');
    hy = ylabel('Latitude');
    set(gca,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    print('figPinfishSST_hianom',['-d' figfmt],'-r300');
    close(fhianom)

    floanom = figure('unit','inches');
    pcolor(lonMAB,latMAB,sstMABlo-sstMAB); shading flat;
    set(gca,'box','on','tickdir','out')
    caxis([-3 3]); colorbar;
    colormap(jet(6/0.1));
    ht = title('SST anomaly of low Pinfish years');
    hx = xlabel('Longitude');
    hy = ylabel('Latitude');
    set(gca,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    print('figPinfishSST_loanom',['-d' figfmt],'-r300');
    close(floanom)
end %if

%hilo
if ismember(plotflag, [4 5 6 7])
    fhilo = figure('unit','inches');
    set(fhilo,'paperposition',[0 0 10 8]);
    hp = pcolor(lonMAB,latMAB,sstMABhi-sstMABlo); shading flat;
    set(gca,'box','on','tickdir','out');
    caxis([-3 3]);
    colorbar;
    colormap(jet(6/0.1));
    ht = title('SST difference between high and low Pinfish years');
    hx = xlabel('Longitude');
    hy = ylabel('Latitude');
    set(gca,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    print('figPinfishSST_hilo',['-d' figfmt],'-r300');
    close(fhilo)
end %if

%yearlyanom
if ismember(plotflag, [8 9 10 12])
    fSSTmap1 = figure('unit','inches')
    set(fSSTmap1,'paperposition',[0 0 10 8]);
    for ii = 1:size(sstMAByr,3)
        hp = pcolor(lonMAB,latMAB,sstMAByr(:,:,ii)-sstMAB); shading flat;
        set(gca,'box','on','tickdir','out')
        caxis([-4 4]); colorbar;
        colormap(jet(6/0.1));
        ht = title(['SST anomaly for ' num2str(tyr(ii))]);
        hx = xlabel('Longitude');
        hy = ylabel('Latitude');
        set(gca,'fontsize',18,'fontweight','bold');
        set(ht,'fontsize',20,'fontweight','bold');
        set(hx,'fontsize',18,'fontweight','bold');
        set(hy,'fontsize',18,'fontweight','bold');
        print(['figPinfishSST_anom' num2str(tyr(ii))],['-d' figfmt],'-r300');
        delete(hp)
    end %for
    close(fSSTmap1)
end %if

%year1anom
if ismember(plotflag, [32 33 34 36 40 48])
    ii = size(sstMAByr,3)

    fSSTmap1 = figure('unit','inches')
    set(fSSTmap1,'paperposition',[0 0 10 8]);
    hp = pcolor(lonMAB,latMAB,sstMAByr(:,:,ii)-sstMAB); shading flat;
    set(gca,'box','on','tickdir','out')
    caxis([-4 4]); colorbar;
    colormap(jet(6/0.1));
    ht = title(['SST anomaly for ' num2str(tyr(ii))]);
    hx = xlabel('Longitude');
    hy = ylabel('Latitude');
    set(gca,'fontsize',18,'fontweight','bold');
    set(ht,'fontsize',20,'fontweight','bold');
    set(hx,'fontsize',18,'fontweight','bold');
    set(hy,'fontsize',18,'fontweight','bold');
    print(['figPinfishSST_anom' num2str(tyr(ii))],['-d' figfmt],'-r300');
    delete(hp)
    close(fSSTmap1)
end %if
