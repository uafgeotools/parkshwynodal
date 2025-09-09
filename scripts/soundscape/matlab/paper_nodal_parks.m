% 
% paper_nodal_parks.m
%
%

close all, clc, clear

ddir = '/home/carltape/PROJECTS/nodal/';

bproposal = false;
bomori = false;
bevents = true;

bprint = false;
bwrite = false;

szone = '5V';
tstart = datenum(2019,2,11);
tend = datenum(2019,3,25);
% denali fault crossing parks highway
loncen = -148.8218;
latcen = 63.4516;
spdy = 86400;

t1 = tstart; t2 = tend;
%t1 = datenum(2019,2,22); t2 = t1+1;

%==========================================================================

if bevents
    % read the whole GCMT catalog -- short version of output
    [otime,tshift,hdur,lat,lon,dep,M,M0,Mw,eid] = readCMT([t1 t2],[],[0 10]);
    display_eq_list([],otime,lon,lat,dep,Mw);
    
%    1      1 otime 2019-02-12 12:34:19.200 lon  146.05 lat   19.01 dep 139.00 km M  6.06
%    2      2 otime 2019-02-14 19:57:12.600 lon  -35.87 lat   35.26 dep  16.80 km M  6.16
%    3      3 otime 2019-02-17 14:35:58.100 lon  152.13 lat   -3.33 dep 372.20 km M  6.38
%    4      4 otime 2019-02-22 10:17:28.000 lon  -77.09 lat   -2.26 dep 121.10 km M  7.49
%    5      5 otime 2019-03-01 08:50:49.000 lon  -70.14 lat  -14.78 dep 273.60 km M  7.04
%    6      6 otime 2019-03-02 03:22:57.900 lon  146.86 lat   41.98 dep  35.20 km M  6.01
%    7      7 otime 2019-03-06 15:46:20.500 lon -177.58 lat  -31.92 dep  12.00 km M  6.40
%    8      8 otime 2019-03-08 15:06:16.400 lon  126.20 lat   10.35 dep  43.30 km M  6.06
%    9      9 otime 2019-03-10 08:12:29.300 lon -178.61 lat  -17.81 dep 580.70 km M  6.27
%   10     10 otime 2019-03-10 12:48:03.300 lon  152.09 lat  -10.13 dep  12.00 km M  6.05
%   11     11 otime 2019-03-15 05:03:57.000 lon  -65.90 lat  -17.74 dep 380.50 km M  6.35
%   12     12 otime 2019-03-20 15:24:02.200 lon  167.48 lat  -15.58 dep 129.70 km M  6.27
%   13     13 otime 2019-03-23 19:21:20.500 lon  -76.31 lat    4.62 dep 128.80 km M  6.08
%   14     14 otime 2019-03-24 04:37:39.100 lon  126.36 lat    1.77 dep  41.90 km M  6.15
    
    % custom file (before AEC Total was updated)
    %ifile = '/home/carltape/dwrite/usgs/query_20190602_F3TN.csv';
    %[otime,lon,lat,dep,mag] = read_eq_usgs(ifile);
    %[otime,lat,lon,dep,mag,isubset] = get_seis_subset(otime,lat,lon,dep,mag,[t1 t2],[],[]);
    [otime,lon,lat,dep,mag,eid,etype] = read_eq_AEC([t1 t2],[],[]);
    
    Mmin = 0;
    Mmin = -10;
    if 1==1
        ax2 = [-151 -148 62 65];
        [otime,lat,lon,dep,mag,isubset] = get_seis_subset(otime,lat,lon,dep,mag,[t1 t2],ax2,[Mmin 10]);
    else
        xedge = [0 360];
        yedge = [-10 300];
        zedge = [0 200 300 1000];
        nkeep = [10 10 10];
        mmin = [3 4 4];
        bdisplayoutput = false;
        [inds,otime,lon,lat,dep,mag] = seis_decluster_circ(xedge,yedge,zedge,...
            loncen,latcen,nkeep,mmin,otime,lon,lat,dep,mag,bdisplayoutput);
    end
    
    isort = [];
    %[~,isort] = sort(mag,'descend');
    display_eq_list(isort,otime,lon,lat,dep,mag);
    
%    1      1 otime 2019-02-18 17:02:46.332 lon -149.93 lat   61.46 dep  40.91 km M  4.36
%    2      2 otime 2019-02-19 06:59:07.339 lon -149.30 lat   61.88 dep   6.91 km M  3.33
%    3      3 otime 2019-02-22 00:28:20.398 lon -149.15 lat   55.63 dep  11.49 km M  4.60
%    4      4 otime 2019-02-23 04:45:11.130 lon -148.56 lat   56.40 dep  11.73 km M  4.40
%    5      5 otime 2019-02-25 18:22:30.737 lon -149.66 lat   62.79 dep  74.48 km M  3.12
%    6      6 otime 2019-03-01 03:41:04.279 lon -148.58 lat   56.37 dep  11.82 km M  4.20
%    7      7 otime 2019-03-03 14:38:12.935 lon -148.69 lat   56.28 dep  10.00 km M  4.20
%    8      8 otime 2019-03-06 00:34:26.853 lon -148.70 lat   62.05 dep  39.64 km M  3.25
%    9      9 otime 2019-03-06 21:33:13.991 lon -157.22 lat   66.31 dep   9.09 km M  5.30
%   10     10 otime 2019-03-07 22:11:51.394 lon -155.07 lat   58.59 dep 130.85 km M  5.00
%   11     11 otime 2019-03-09 22:26:09.184 lon -145.12 lat   63.25 dep   5.59 km M  3.69
%   12     12 otime 2019-03-09 23:39:58.332 lon -147.74 lat   64.55 dep  26.50 km M  3.61
%   13     13 otime 2019-03-10 05:16:05.748 lon -150.52 lat   62.88 dep  89.22 km M  3.08
%   14     14 otime 2019-03-10 22:40:10.677 lon -141.59 lat   59.81 dep   0.02 km M  4.23
%   15     15 otime 2019-03-11 00:07:18.295 lon -141.63 lat   59.79 dep  13.82 km M  5.00
%   16     16 otime 2019-03-12 15:52:51.565 lon -157.21 lat   66.28 dep   8.15 km M  4.49
%   17     17 otime 2019-03-13 02:42:25.318 lon -144.50 lat   69.52 dep  13.00 km M  4.29
%   18     18 otime 2019-03-14 14:16:42.542 lon -149.12 lat   63.82 dep 119.31 km M  3.15
%   19     19 otime 2019-03-15 05:58:36.024 lon -151.13 lat   63.28 dep   9.03 km M  3.93
%   20     20 otime 2019-03-19 02:39:35.423 lon -150.42 lat   63.18 dep 116.34 km M  3.42
%   21     21 otime 2019-03-23 15:14:45.096 lon -149.86 lat   61.53 dep  47.39 km M  4.06
    
    ftag = '/home/carltape/REPOSITORIES/GEOTOOLS/python_util/pysep_dev/input/denali_parks';
    if bwrite, write_seis_obspy(ftag,lon,lat,dep,mag,otime); end
    
    if 0==1
        % find ALL events within 3 minutes of each of these
        for ii=1:length(otime)
            t1 = otime(ii) + 1/spdy;
            t2 = t1 + 180/spdy;
            [otimex,lonx,latx,depx,magx] = read_eq_AEC([t1 t2],[],[]);
            display_eq_list([],otimex,lonx,latx,depx,magx);
        end

        % ALL events in the catalog
        [otime,lon,lat,dep,mag] = read_eq_AEC([t1 t2],[],[]);
        display_eq_list([],otime,lon,lat,dep,mag);
    end
    
end
    
%==========================================================================

if bproposal
    addpath(ddir);  % for kml2struct.m
    stag = 'dead_mouse'; N = 94;
    %stag = 'parks_hwy'; N = 306;
    %stag = 'parks_hwy_epi';
    kmlfile = [ddir stag '.kml'];   % digitized in google earth
    
    k = kml2struct(kmlfile);
    klon = k.Lon;
    klat = k.Lat;
    % order points from south to north
    if max(klat)==klat(1), klat=flipud(klat); klon=flipud(klon); end
    
    % utm
    [kx,ky] = utm2ll(klon,klat,szone,0);
    %kx = klon; ky = klat;
    
    % interpolate the digitized curve (probably better done in UTM)
    smethod = 'linear';
    %smethod = 'spline';
    t = linspace(0,1,N);
    [pt,dudt,fofthandle] = interparc(t,kx,ky,smethod);
    kxi = pt(:,1);
    kyi = pt(:,2);
    plot(kxi,kyi,'r');
    
    % check distance between points
    rdist = sqrt( (kxi(1:N-1) - kxi(2:N)).^2 + (kyi(1:N-1) - kyi(2:N)).^2 );
    figure; plot_histo(rdist,[700:50:1300]);
    
    [kloni,klati] = utm2ll(kxi,kyi,szone,1);
    
    figure; hold on;
    plot(klon,klat,'bo-','markersize',10,'markerfacecolor','b');
    plot(kloni,klati,'ro','markersize',10);
    
    if bwrite
        write_xy_points([ddir stag],klon,klat,'%12.6f %12.6f');
        write_xy_points([ddir stag '_' smethod],kloni,klati,'%12.6f %12.6f');
        kmlwrite([ddir stag '_' smethod],klati,kloni,'Name',numvec2cell(1:N));
    end
    
%     ifile3 = '/home/admin/share/datalib/cities_roads_etc/roads/alaska_roads.xy';
%     [data3,ilabs,labels] = read_delimited_file(ifile3,2);
%     rlon = data3(:,1);
%     rlat = data3(:,2);
%     figure; hold on;
%     plot(rlon,rlat,'k','linewidth',2);
%     
%     ax = [-150.2599 -148.7105   61.3284   64.9499];
%     axis(ax);
%     ilab = numvec2cell([1:length(rlon)]);
%     text(rlon,rlat,ilab);
%    
%     lonS = -149.9911; latS = 61.7697;   % Willow
%     lonN = -149.0931; latN = 64.5639;   % Nenana
%     
%     [~,iN] = min(distance(latN*ones(size(rlat)),lonN*ones(size(rlat)),rlat,rlon));
%     [~,iS] = min(distance(latS*ones(size(rlat)),lonS*ones(size(rlat)),rlat,rlon));
%     plot(rlon(iN),rlat(iN),'rp');
%     plot(rlon(iS),rlat(iS),'rp');
end

if 0==1
    idir = '/home/carltape/PROJECTS/earthquakes/20181130_anchorage/';
    ifile = [idir 'CMTSOLUTION_201811301729A'];
    [otime,tshift,hdur,lat,lon,dep,M,eid,elabel] = read_CMTSOLUTION(ifile);
    slabel = 'anchorage2018';
    write_psmeca([idir slabel],otime,lat,lon,dep,M,eid);
end

if bomori
    % aftershocks for 2018 Anchorage earthquake
    MMIN = 2;
    ifile = '/home/carltape/dwrite/usgs/query_anchorage2018_v2.csv';
    stag = sprintf('M >= %.1f, 2018',MMIN);
    [otime,lon,lat,dep,mag] = read_eq_usgs(ifile);
    [~,isort] = sort(mag,'ascend');
    display_eq_summary(otime,lon,lat,dep,mag);
    display_eq_list(isort,otime,lon,lat,dep,mag);
    
    error
    
    DREF = datenum(2018,11,1);
    DQ = datenum(2018,11,30,17,29,29);
    DMAX = datenum(2019,4,1);
    doy = floor(otime - DREF);
    inds = find(mag >= MMIN);
    figure; scatter(lon(inds),lat(inds),4^2,doy(inds),'filled');
    
    figure;
    for kk=1:2
        subplot(4,1,kk); hold on;
        [N,Nplot,centers] = plot_histo(doy(inds),[DREF:DMAX]-DREF,1);
        %title(stag);
        xlabel(sprintf('Day since %s',datestr(DREF,29)));

        % omori curve
        N0 = 0;
        K = -1;
        NDAY = round(DMAX-DREF);
        t = [-K:1:NDAY];
        C = 300;
        p = 0.9;
        n = C./(K+t).^p + N0;
        tplot = t + (DQ-DREF);
        plot(tplot,n,'r--','linewidth',2);
        xvec = [datenum(2019,2,10) datenum(2019,3,27)] - DREF;
        [xmat,ymat] = vertlines(xvec,0,400);
        plot(xmat,ymat,'k');
        if kk==1, ylim([0 400]); else ylim([0 10]); end
    end
    
    if bprint, orient tall; print(gcf,'-dpsc',sprintf('%somori',ddir)); end
    
    % frequency-magnitude plots
    dmag = 0.1;
    inds = find(otime > DQ);    % before mainshock
    tran_yr = otime(inds);
    [N,Ninc,Medges] = seis2GR(mag(inds),dmag);
    [b1,a1,a1i] = GR2plot(N,Ninc,Medges,tran_yr,[2 3]);
    inds = find(otime < DQ);    % after mainshock
    tran_yr = otime(inds);
    [N,Ninc,Medges] = seis2GR(mag(inds),dmag);
    [b1,a1,a1i] = GR2plot(N,Ninc,Medges,tran_yr,[2 3]);
end
    
%==========================================================================

if 0==1
    idir = '/home/carltape/PROJECTS/earthquakes/20181130_anchorage/';
    ifile = [idir 'CMTSOLUTION_201811301729A'];
    [otime,tshift,hdur,lat,lon,dep,M,eid,elabel] = read_CMTSOLUTION(ifile);
    slabel = 'anchorage2018';
    write_psmeca([idir slabel],otime,lat,lon,dep,M,eid);
end

%==========================================================================
