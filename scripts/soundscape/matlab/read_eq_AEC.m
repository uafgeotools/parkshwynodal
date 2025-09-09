function [otime,lon,lat,dep,mag,eid,etype] = read_eq_AEC(oran,ax3,Mwran)
%READ_EQ_AEC Read seismicity catalog from the Alaska Earthquake Center.
%
% Option arguments to extract a subset based on origin time, region, and
% magnitude.
%
% WARNING 1: THE TOTAL CATALOG CONTAINS ONLY OFFICIAL LOCATIONS AND IS
%            MANY MONTHS BEHIND REAL TIME.
% WARNING 2: THE MATLAB VERSION, LOADED HERE, NEEDS TO BE RE-SAVED EACH
%            TIME THAT TOTAL IS UPDATED.
%
% To check when Total was updated:
%   ls -ltr /aec/db/catalogs/final/Total*
% BETTER (check date of event at the bottom):
%   dbe /aec/db/catalogs/final/Total.origin
% See instructions to update the Total.mat file here:
%   /home/admin/share/datalib/seismicity/regional/alaska/AEC/README
%
% To write subsets of this catalog for GMT plotting, see write_cat_subsets.m
%
% calls get_seis_subset.m
%

% open a custom file that includes events more recent than those in Total
% see /home/admin/share/datalib/seismicity/regional/alaska/AECex/
buse_extended = false;

% custom extended file
if buse_extended
    idir = '/home/admin/share/datalib/seismicity/regional/alaska/AECex/';
    cfile = strcat(idir,'Total_merged');
    ctag = ' and extended';
else
    idir = '/home/admin/share/datalib/seismicity/regional/alaska/AEC/';
    cfile = '/aec/db/catalogs/final/Total';
    ctag = '';
end

ftag = 'Total_read_eq_AEC';
if nargout==7
    ftag = 'Total_read_eq_AEC_etype';
end
fullfile = [idir ftag '.mat'];
eid = [];

if ~or(nargin==3,nargin==0), error('input arguments must be 0 or 3'); end

if exist(fullfile,'file')
    disp('loading the full AEC catalog (matlab-saved Total), minus glacier (E) and quarry (Q) events');
    load(fullfile);
else
    disp(sprintf('reading the full%s AEC catalog (Total)',ctag));

    % UAF/GI black box ('help aeic_catalog')
    %/home/admin/share/matlab/PACKAGES/GISMO/uaf_internal/AEIC_AVO/+aeic_catalog/get_total.m
    catalog = aeic_catalog.get_total(cfile);

    lat = catalog.lat;
    lon = catalog.lon;
    dep = catalog.depth;
    otime = catalog.time;
    mag = catalog.prefMagnitude;
    n = length(lon);
    
    % event types
    if 0==1
        etype = catalog.etype;
        uetype = unique(etype);
        for kk=1:length(uetype)
            imatch = find(strcmp(etype,uetype{kk}));
            stit = sprintf('%7i / %7i = %.3f are event type %s',length(imatch),n,length(imatch)/n,uetype{kk});
            disp(stit);
            figure; ax0 = [190 230 56 72];
            alaska_basemap(ax0,0,0,0)
            plot(wrapTo360(lon(imatch)),lat(imatch),'.');
            title(stit); axis(ax0);
        end
    end
    
    % exclude glacier events and quarry events
    % (note that unknown X events are left in)
    etype = catalog.etype;
    ecut = {'G','Q'};
    icut = [];
    for kk=1:length(ecut)
        imatch = find(strcmp(etype,ecut{kk}));
        icut = [icut ; imatch];
        disp(sprintf('%i events of type %s to cut',length(imatch),ecut{kk}));
    end
    lat(icut) = [];
    lon(icut) = [];
    dep(icut) = [];
    otime(icut) = [];
    mag(icut) = [];
    etype(icut) = [];
    
    % Saving the etype makes it slower to load (2.0 s vs 0.2 s);
    % in most cases the user does not need to know the etype.
    disp('saving the AEC catalog (matlab-saved Total), minus glacier (E) and quarry (Q) events');
    sfile = sprintf('./%s_%s',ftag,datestr(now,'yyyymmdd'));
    if nargout==7
        save(sfile,'otime','lon','lat','dep','mag','etype');
    else
        save(sfile,'otime','lon','lat','dep','mag');
    end
end

% optional: extract a subset
if nargin == 3
    [otime,lat,lon,dep,mag,isubset] = get_seis_subset(otime,lat,lon,dep,mag,oran,ax3,Mwran);
    
    % subset the extra variables
    if isempty(isubset)
        eid = [];
        return
    end
    eid = otime2eid(otime);
    if nargout==7, etype = etype(isubset); end
end

%==========================================================================
% EXAMPLES

if 0==1
    % load full catalog (no event ID)
    % (it takes much longer to load the etype)
    tic, [otime,lon,lat,dep,mag] = read_eq_AEC; toc
    tic, [otime,lon,lat,dep,mag,eid,etype] = read_eq_AEC; toc
    display_eq_summary(otime,lon,lat,dep,mag);
    
    %--------------
    % exact origin time for a target event
    [otime,lon,lat,dep,mag] = read_eq_AEC([datenum(2016,5,18,3,25,48) 10],[],[]);
    sotime = datestr(otime,'yyyy-mm-dd HH:MM:SS.FFF')
    otime2eid(otime)
    
    %--------------
    load coast;
    clat = lat; clon = long;
    
    %% GCMT largest (intraplate) earthquakes since 2002 Denali (SSA 2022 talk)
    ax2 = [-180 -141 56 75];
    oran = [datenum(2002,11,2) today];
    [otime,tshift,hdur,lat,lon,dep,M,M0,mag,eid] = readCMT(oran,ax2,[5.95 10]);
    display_eq_summary(otime,lon,lat,dep,mag);
    [~,isort] = sort(otime);
    display_eq_list(isort,otime,lon,lat,dep,mag);

    % Labeling [XX/YY]:
    % XX is the plate whose boundary it is inside
    % YY is the plate it is within
    %    1      1 otime 2002-11-03 22:13:28.000 lon -144.89 lat   63.23 dep  15.00 km M  7.849 -- Denali
    %    2      2 otime 2009-03-30 07:13:08.800 lon -152.64 lat   56.26 dep  12.00 km M  6.027 -- megathrust
    %    3      3 otime 2010-04-30 23:11:49.000 lon -177.71 lat   60.45 dep  15.00 km M  6.483 -- [NA/NA] Bering
    %    4      4 otime 2010-04-30 23:16:32.600 lon -177.47 lat   60.29 dep  16.00 km M  6.352 -- [NA/NA] Bering
    %    5      5 otime 2012-11-12 20:42:16.900 lon -142.97 lat   57.70 dep  16.00 km M  6.344 -- [PA/PA] Transition fault
    %    6      6 otime 2014-09-25 17:51:22.400 lon -151.78 lat   62.00 dep 109.20 km M  6.287 -- [NA/PA] intraslab ~Susitna
    %    7      7 otime 2015-05-29 07:00:13.800 lon -156.52 lat   56.71 dep  81.10 km M  6.798 -- [PA/PA] intraslab ~Kodiak
    %    8      8 otime 2015-07-29 02:36:01.600 lon -153.15 lat   60.02 dep 129.80 km M  6.366 -- [PA/PA] intraslab ~Iniskin
    %    9      9 otime 2016-01-24 10:30:37.400 lon -153.27 lat   59.75 dep 110.70 km M  7.115 -- [NA/PA] Iniskin
    %   10     10 otime 2016-04-02 05:50:08.000 lon -157.66 lat   57.22 dep  12.00 km M  6.221 -- [NA/PA] intraslab ~Kodiak
    %   11     14 otime 2018-01-23 09:32:04.000 lon -149.12 lat   56.22 dep  33.60 km M  7.928 -- [PA/PA] offshore Kodiak**
    %   12     11 otime 2018-08-12 14:59:02.200 lon -144.78 lat   69.74 dep  12.00 km M  6.430 -- [NA/NA] Kaktovik
    %   13     12 otime 2018-08-12 21:15:04.100 lon -144.23 lat   69.68 dep  15.30 km M  6.015 -- [NA/NA] Kaktovik
    %   14     16 otime 2018-11-30 17:29:34.700 lon -150.02 lat   61.49 dep  48.20 km M  7.054 -- [NA/PA] intraslab ~Nelchina**
    %   15     15 otime 2021-05-31 06:59:57.600 lon -148.41 lat   62.57 dep  66.30 km M  6.098 -- [NA/PA] intraslab ~Nelchina
    %   16     17 otime 2021-10-11 09:10:30.000 lon -156.63 lat   56.45 dep  66.60 km M  6.908 -- [PA/PA] intraslab ~Kodiak
    %   17     13 otime 2021-12-21 22:42:16.200 lon -153.41 lat   60.16 dep 154.10 km M  5.965 -- [NA/PA] intraslab ~Iniskin

    figure; hold on; plot(clon,clat);
    plot(lon,lat,'.'); stlab = [];
    for ii=1:length(lon), stlab{ii} = sprintf('%s M%.1f',datestr(otime(ii),'yyyymmdd'),mag(ii)); end
    text(lon,lat,stlab); axis(ax2);
    
    %--------------
    %%
    
    ax3 = [-151 -145   63.5 66 -10 700];      % interior
    % ax3 = [-151 -147.5 63.5 65.5 -10 700];    % MFSZ
    % ax3 = [-150.5 -148.2 63.8 65.3 -10 700];  % MFSZ (zoomed)
    % ax3 = [-150.5 -149.7 63.8 64.2 -10 700];  % MFSZ (south)
    % ax3 = [-149.7 -148.2 64.4 65.3 -10 700];  % MFSZ (north,southwest,east,west)
    oran = [datenum(2000,1,1) datenum(2100,1,1)];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag,eid] = read_eq_AEC(oran,ax3,Mwran);
    n = length(otime);
    
    figure; scatter(lon,lat,4^2,year(otime),'filled');
    axis tight; colorbar; caxis([1990 2010]); 
    title('Colored by year');
    figure; scatter(lon,lat,4^2,dep,'filled');
    axis tight; colorbar; caxis([0 25]);
    title('Colored by depth, km');
    
    %--------------
    % southern Alaska 
    ax3 = [-158 -142 55 63 -10 700];
    oran = [datenum(1990,1,1) datenum(2100,1,1)];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag,eid] = read_eq_AEC(oran,ax3,Mwran);
    n = length(otime);
    
    figure; scatter(lon,lat,4^2,dep,'filled');
    axis tight; colorbar; caxis([0 25]);
    title('Colored by depth, km');
    figure; scatter3(lon,lat,dep,4^2,dep,'filled');colorbar;
    xlabel('lon'),    ylabel('lat'),    zlabel('depth');
    
    %----------------------------------
    % 1995 MFSZ aftershocks
    ax3 = [-151 -147.5 63.5 65.5 -10 700];
    oran = [datenum(1995,10,6) datenum(1996,10,6)];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    n = length(otime);
    figure; scatter(lon,lat,4^2,dep,'filled');
    axis tight; colorbar; caxis([0 25]);
    title('Colored by depth, km');
    imain = find(mag==6);
    lon(imain), lat(imain), dep(imain)
    
    %--------------
    
    % intraplate deep crustal seismicity
    dran = [15 35]; Mwran = [2 10]; yrmin = 2000;
    %dran = [20 35]; Mwran = [0 10]; yrmin = 2000;
    %dran = [-10 10]; Mwran = [0 10]; yrmin = 2003;
    ax3 = [-170 -140 62 70 dran];
    oran = [datenum(yrmin,1,1) now];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    n = length(otime);
    figure; alaska_basemap(ax3(1:4),0,0,0);
    %hold on; plot(clon,clat); axis(ax3(1:4));
    scatter(wrapTo360(lon),lat,4^2,dep,'filled');
    colorbar; caxis(dran);
    title('Colored by depth, km');
    
    %--------------
    
    % clusters of deep seismicity
    dran = [-10 700];
    %dran = [20 35];
    ax3 = [201 204 66 66.7 dran];
    ax3 = [216 219 66 66.6 dran];
    ax3 = [207 214 65.3 66.7 dran];
    oran = [datenum(1990,1,1) now];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    figure; scatter(lon,lat,6^2,dep,'filled');  
    title('colored by depth, km');
    caxis([0 30]);
    colorbar; axis tight;
    
    %--------------
    
    % Fairbanks seismicity -- time migration?
    ax3 = [-148 -147 64.7 65.1 -10 700];
    oran = [datenum(1990,1,1) datenum(2200,1,1)];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    figure; scatter(lon,lat,6^2,year(otime),'filled');  
    title('Fairbanks region seismicity, colored by year');
    colorbar; axis tight;
    
    %--------------
    
    % Cook Inlet crustal seismicity
    ax3 = [-153.5 -149 59 61.75 -10 20];
    oran = [datenum(1990,1,1) datenum(2011,1,1)];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    figure; scatter(lon,lat,6^2,dep,'filled');  
    title('Cook Inlet crustal seismicity, colored by depth');
    colorbar; axis tight;
    
    %--------------
    
    % Cook Inlet seismicity
    ax3 = [-153.5 -149 59 61 -10 30];
    oran = [datenum(1990,1,1) datenum(2011,1,1)];
    Mwran = [2 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    figure; scatter(lon,lat,6^2,dep,'filled');  
    title('Cook Inlet crustal seismicity, colored by depth');
    colorbar; axis tight;
    
    %--------------
    
    % Seismicity around Skagway eq
    ax3 = [-139 -134 58 61 -10 40];
    oran = [datenum(1900,1,1) today];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    figure; scatter(wrapTo360(lon),lat,6^2,dep,'filled');  
    title('Skagway seismicity 1990-2017, colored by depth');
    colorbar; axis tight;
    hold on
    alaska_basemap(ax3);
    plot(wrapTo360(-136.7163),59.8184,'p','MarkerSize',20,'MarkerFaceColor','c','MarkerEdgeColor','k')
    plot(wrapTo360(-136.6618),59.8522,'p','MarkerSize',20,'MarkerFaceColor','c','MarkerEdgeColor','k')
    %--------------
    
    % Cook Inlet seismicity -- day-time vs night-time
    % Alaska time is GMT - 9
    ax3 = [-153.5 -149 59 61.75 -10 20];
    oran = [datenum(1990,1,1) datenum(2011,1,1)];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    n = length(otime);
    inight = find(and(hour(otime) > 9, hour(otime) < 21));
    iday = setdiff([1:n]',inight);
    
    figure; nr=2; nc=1; 
    subplot(nr,nc,1); plot(lon(inight),lat(inight),'.');
    title(sprintf('%i / %i events during the 12-hour night',length(inight),n));
    subplot(nr,nc,2); plot(lon(iday),lat(iday),'.');
    title(sprintf('%i / %i events during the 12-hour day',length(iday),n));
    
    %--------------
    
    % Cook Inlet seismicity
    ax3 = [-152 -150 61 62 -10 20];
    oran = [datenum(1990,1,1) datenum(2011,1,1)];
    Mwran = [0 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    figure; scatter(lon,lat,6^2,dep,'filled');  
    title('Cook Inlet crustal seismicity, colored by depth');
    colorbar; axis tight;
    
    % Talachulitna seismicity (see Flores and Doser, 2005)
    ax3 = [-151.5 -151.1 61.25 61.90 -10 40];
    oran = [datenum(1900,1,1) datenum(2018,1,1)];
    Mwran = [1.5 10];
    Mwran = [];
    %[otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    [otime,lon,lat,dep,M,M0,mag,eid] = read_mech_AECfp(oran,ax3,Mwran);
    figure; scatter(lon,lat,6^2,dep,'filled'); axis(ax3(1:4)); 
    title('Talachulitna crustal seismicity, colored by depth');
    colorbar; axis tight;
    [~,isort] = sort(mag,'descend')
    display_eq_list(isort,otime,lon,lat,dep,mag);
    for ii=1:length(otime), figure; plot_beachballs(M(:,ii)); end
    PAall = CMT2pa(M)
    PAall(:,[6 5 2 1])
    
    %% Castle Mountain 1984
    otar = datenum(1984,8,14,1,2,0);
    [otime,lon,lat,dep,mag] = read_eq_AEC([otar 60],[],[]);
    display_eq_list([],otime,lon,lat,dep,mag);
    % Sutton 1996
    otar = datenum(1996,11,11,10,52,0);
    [otime,lon,lat,dep,mag] = read_eq_AEC([otar 60],[],[]);
    display_eq_list([],otime,lon,lat,dep,mag);
    [otime,lon,lat,dep,M,M0,mag,eid] = read_mech_AECfp([otar 60],[],[]);
    figure; plot_beachballs(M);
    %% crustal near Anchorage 1997
    otar = datenum(1997,5,6,1,31,0);
    [otime,lon,lat,dep,mag] = read_eq_AEC([otar 60],[],[]);
    display_eq_list([],otime,lon,lat,dep,mag);
    [otime,lon,lat,dep,M,M0,mag,eid] = read_mech_AECfp([otar 60],[],[]);
    figure; plot_beachballs(M);
    [otime,tshift,hdur,lat,lon,dep,M,M0,Mw] = readCMT('B050697B');
    figure; plot_beachballs(M);
    
    %--------------
    
    %% slab seismicity in Alaska (depth extent of slab is digitized into ddir)
    dmin = 40;  % 60
    Mmin = 1;  % 2
    %ax3 = [-165 -140 49 70 dmin 700];
    ax3 = [160 220 49 70 dmin 700];   % western part
    oran = [datenum(1990,1,1) datenum(2013,1,1)];
    Mwran = [Mmin 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    [~,isort] = sort(dep,'descend');
    % plot figure for digitizing
    figure; scatter(wrapTo360(lon(isort)),lat(isort),4^2,dep(isort),'filled');  
    title('Alaska slab (and lower crustal) seismicity, colored by depth');
    caxis([dmin 140]); colorbar; axis tight; grid on;
    xlabel('Longitude'); ylabel('Latitude');
    % plot slab extent boundaries
    ddir = '/home/admin/share/datalib/seismicity/regional/alaska/slabextent/';
    d = load([ddir 'slab_extent.lonlat']);
    hold on; plot(wrapTo360(d(:,1)),d(:,2),'k--','linewidth',2);
    d = load([ddir 'slab_extent_wrangell.lonlat']);
    plot(wrapTo360(d(:,1)),d(:,2),'k--','linewidth',2);
    % plot yakutat boundary
    %d = load(['/home/carltape/gmt/misc/YAK_extent.lonlat']);
    %plot(wrapTo360(d(:,1)),d(:,2),'k','linewidth',2);
    
    %--------------
    
    % Wrangell subduction zone (see Wang and Tape supplement)
    dmin = 40; dmax = 120; Mmin = 1;
    wrangell_mean_lat = 62.3;
    Mwran = [Mmin 10];
    ax3 = [214.5 219 59 63.5 -10 700];
    oran = [datenum(1990,1,1) datenum(2013,1,1)];
    %oran = [datenum(1990,1,1) now];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    islab = find(dep >= dmin);
    % plot figure
    figure; nr=2; nc=1;
    % slab events in map view
    subplot(nr,nc,1);
    % slab extent boundaries
    ddir = '/home/admin/share/datalib/seismicity/regional/alaska/slabextent/';
    d = load([ddir 'slab_extent_wrangell.lonlat']);
    plot(wrapTo360(d(:,1)),d(:,2),'r','linewidth',4);
    % basemap
    alaska_basemap(ax3(1:4),0,0,0);
    plot(ax3(1:2),wrangell_mean_lat*[1 1],'r--');
    %hold on; plot(clon,clat); axis(ax3(1:4));
    scatter(wrapTo360(lon(islab)),lat(islab),8^2,dep(islab),'filled');
    scatter(wrapTo360(lon(islab)),lat(islab),8^2,'ko');
    colorbar; caxis([dmin dmax]); axis equal, axis(ax3(1:4));
    title(sprintf('(a) %i events below %.0f km, colored by depth',length(islab),dmin));
    % seismicity in vertical cross section
    subplot(nr,nc,2); hold on;
    set(gca,'ydir','reverse');
    plot(lat,dep,'.');
    scatter(lat(islab),dep(islab),8^2,dep(islab),'filled');
    scatter(lat(islab),dep(islab),8^2,'ko');
    caxis([dmin dmax]);
    ylim([0 dmax]);
    plot(ax3(3:4),[dmin dmin],'r--');
    plot(wrangell_mean_lat*[1 1],[0 dmax],'r--');
    [xf,yf,mf,bf,rms] = linefit([61.8 62.6],[dmin dmax]);
    plot(xf,yf,'r--');
    xlabel('Latitude'); ylabel('Depth, km');
    title(sprintf('(b) AEC catalog, %.0f-%.0f, M >= %.1f: all %i events in map region of (a)',...
        year(oran(1)),year(oran(2)),Mmin,length(lat)));

    %--------------
    
    % ARCTIC events
    ax3 = [-156 -144 64 72 -10 700];
    oran = [datenum(2005,5,1) datenum(2007,7,1)];
    Mwran = [4 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);
    n = length(lon);
    [rlon,rlat,relev,rburial,stnm,netwk] = read_station_ARCTIC;
    nrec = length(rlon);
    hold on; scatter(lon,lat,6^2,dep,'filled'); grid on;
    title('ARCTIC seismicity, colored by depth');
    caxis([0 30]); colorbar; axis tight;
    %[~,isort] = sort(mag,'descend');
    [~,isort] = sort(lat,'ascend');
    display_eq_list(isort,otime,lon,lat,dep,mag);
    % closest station to each event
    for ii=1:n
        [dmin,imin] = min(deg2km(distance(lat(ii)*ones(nrec,1),lon(ii)*ones(nrec,1),rlat,rlon)));
        disp(sprintf('%s is %5.1f km from %s event at (%6.2f, %6.2f)',...
            stnm{imin},dmin,datestr(otime(ii),29),lon(ii),lat(ii)));
    end
    
    %--------------  
    
    % events in the best-coverage MOOS region
    ax0 = [-148.4048 -146.7565   64.2776   65.4013];
    oran = [datenum(2009,5,1) datenum(2011,1,1)];
    Mwran = [3 10];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,[ax0 -10 700],Mwran);
    display_eq_list(1:length(otime),otime,lon,lat,dep,mag);
    
    %--------------  
    
    % Cook Inlet crustal faults, MOOS time period
    oran = [datenum(2007,8,15) datenum(2009,8,15)];
    ax0 = [-152.5 -149.5 60 61.5];
    ax3 = [ax0 -10 30]; 
    Mwran = [2.8 10];  % try M > 2
    [otime,slat,slon,sdep,mag,eid] = read_eq_AEC(oran,ax3,Mwran);
    %[~,isort] = sort(sdep,'ascend');
    [~,isort] = sort(mag,'descend');
    display_eq_list(isort,otime,slon,slat,sdep,mag,eid);
    
    % same events, AEIC first motion solution
    [otime,slat,slon,sdep,M,M0,mag,eid] = read_mech_AEICfp(oran,ax3,Mwran);
    [~,isort] = sort(mag,'descend');
    display_eq_list(isort,otime,slon,slat,sdep,mag,eid);
    
    % same events, AEIC moment tensor
    [otime,slat,slon,sdep,M,M0,mag,eid,depc] = read_mech_AEIC(oran,ax3,Mwran);
    [~,isort] = sort(mag,'descend');
    display_eq_list(isort,otime,slon,slat,sdep,mag,eid,depc);
    
    % all events in the region, all time
    oran = datenum(2000,1,1);
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax0,[]);
    figure; plot_histo(dep,[-10:2:160]); xlabel('Depth, km');
    
    %--------------------
    
    % SCAK region for CAP (Silwal and Tape, 2016 JGR)
    ax0 = [-154 -146 58 62.5];    
    ax3 = [ax0 -10 200];
    oran = [datenum(2007,8,15) datenum(2009,8,15)];  % MOOS
    oran = [datenum(2015,6,1) datenum(2016,6,1)];    % SALMON year 1
    %oran = [datenum(2016,1,24) datenum(2016,6,1)];   % Iniskin aftershocks
    Mwran = [3.5 10];   % SilwalTape2016
    Mwran = [4 10];
    [otime,lon,lat,dep,mag,eid] = read_eq_AEC(oran,ax3,Mwran);
    figure; scatter(lon,lat,6^2,dep,'filled');  
    %title('Events for MT inversion using CAP');
    colorbar; axis tight;
    [~,isort] = sort(mag,'descend');
    display_eq_list(isort,otime,lon,lat,dep,mag);
    figure; plot(otime,mag,'ko','markersize',14,'markerfacecolor','r');
    datetick('x','mm/yy'); xlabel('Date'); ylabel('Magnitude');
    
    %--------------------
    % Iniskin aftershocks
    axiniskin = [-153.5 -152.7 59.5 60];
    ax3 = [axiniskin 40 700];
    oday = datenum(2016,1,24);
    sdur = 250;
    for kk=1:2
        if kk==1
            oran = [oday oday+sdur];       % Iniskin aftershocks 
        else
            oran = [oday oday+sdur]-365;   % background seismicity
        end
        disp(sprintf('%s to %s',datestr(oran(1)),datestr(oran(2))));
        Mwran = [2 10];
        [otime,lon,lat,dep,mag,eid] = read_eq_AEC(oran,ax3,Mwran);
        figure; hold on; scatter(lon,lat,6^2,dep,'filled');  
        [~,imax] = max(mag); plot(lon(imax),lat(imax),'kp','markersize',12);
        title('region of M7.1 Iniskin earthquake');
        colorbar; axis tight;
        [~,isort] = sort(mag,'descend');
        display_eq_list(isort,otime,lon,lat,dep,mag);
        figure; plot(otime,mag,'ko','markersize',14,'markerfacecolor','r');
        datetick('x','mm/yy'); xlabel('Date'); ylabel('Magnitude');
        tplot = otime-oran(1);
        tedge = [0:1:sdur];
        N = plot_histo(tplot,tedge,1,false);
        if kk==1
            figure(10); grid on; plot(tedge(2:end),cumsum(N),'b.-');
        else
            yshift = 460;
            figure(10); hold on; plot(tedge(2:end),cumsum(N)+yshift,'r.-');
            legend('2016','2015','location','northwest'); grid on; 
            plot([0 sdur],yshift*[1 1],'r--');
            xlabel('day since Jan 24');
            ylabel('cumulative number of earthquakes since Jan 24');
        end
    end

    %--------------------
    % northwest Alaska
    Mwran = [4 10];
    %ax3 = [213 221 59 63.5 dmin 700];
    ax3 = [-168 -158 65 75 -10 700];
    oran = [datenum(1960,1,1) now];
    [otime,lon,lat,dep,mag] = read_eq_AEC(oran,ax3,Mwran);

    %--------------------
    
    % for more RECENT EVENTS that are not within Total
    dbName = '/Seis/processing/analyzed/2012_05/analyzed_2012_05_16';
    db = dbopen(dbName,'r');
    db  = dblookup(db,'','origin','','');
    [otime,elat,elon,edep,eml] = ...
         dbgetv(db,'origin.time','origin.lat','origin.lon','origin.depth','origin.ml');
    dbclose(db);
    otime = epoch2datenum(otime);
    % subset for cook inlet
    %axbox = [-156 -146 58 62];
    axbox = [-153.5 -149 59.6 61.8];
    [otime,elat,elon,edep,eml] = get_seis_subset(otime,elat,elon,edep,eml,[],[axbox -10 700],[]);
    %[~,isort] = sort(eml,'descend');
    [~,isort] = sort(otime,'ascend');
    display_eq_list(isort,otime,elon,elat,edep,eml);

    %--------------------

    % compare GCMT catalog magnitudes vs Total
    % pick limits for the catalog searches
    % (ideally you might want to use the entire Total catalog, without any subset)
    oran = [datenum(1976,1,1) now];
    ax2 = [160 240 45 74];
    Mran = [4 10];
    % GCMT (set A)
    [otime1,~,~,lat1,lon1,dep1,M,M0,mag1,eid1] = readCMT(oran,ax2,Mran);
    figure; plot(wrapTo360(lon1),lat1,'.');
    % Total (set B)
    [otime2,lon2,lat2,dep2,mag2,eid2] = read_eq_AEC(oran,ax2,Mran);
    figure; plot(wrapTo360(lon2),lat2,'.');
    % match to Total
    otime_diff_sec_max = 180;
    hdist_diff_km = 200;
    args = [otime_diff_sec_max hdist_diff_km];
    idisplay = 1;
    idisplay2 = 0;
    [imatch1,imatch2,imatch_otime,otime_diff_sec_all,...
        hdist_diff_km_all,vdist_diff_km_all,mag_diff_all,i1_1notin2,i2_2notin1] ...
        = seismicity_match(args,otime1,otime2,...
         lon1,lon2,lat1,lat2,dep1,dep2,mag1,mag2,eid1,eid2,idisplay,idisplay2);
    figure; plot_histo(mag_diff_all,[-1:0.1:3],1); 
    xlabel('M-GCMT minus M-Total');

    % example event that has a big difference between the magnitudes
    %1535 otime1 2013-09-12 15:25:23, Dotime =  3.1 s, Depi = 15.6 km, Ddep = -17.5 km, DMw =  0.95 C201309121525A 20130912152520020
    %1535 2067 otime 2013-09-12 (2013255) 15:25:23 lon -171.15 lat  51.52 dep  12.00 km Mw 5.09 C201309121525A 20130912152520020
    igcmt = 2067;
    otimetar = otime1(igcmt);
    datestr(otimetar,31)
    [otimex,lonx,latx,depx,magx,eidx] = read_eq_AEC([otimetar 10],[],[]);
    display_eq_list([],otimex,lonx,latx,depx,magx);
    disp('GCMT mag:'); mag1(igcmt)
    disp('Total mag:'); magx
    
    %--------------------
    
    %% read the full catalog, including all event types
    idir = '/home/admin/share/datalib/seismicity/regional/alaska/AEC/';
    cfile = '/aec/db/catalogs/final/Total';
    catalog = aeic_catalog.get_total(cfile);
    lat = catalog.lat;
    lon = catalog.lon;
    dep = catalog.depth;
    otime = catalog.time;
    mag = catalog.prefMagnitude;
    etype = catalog.etype;
    n = length(lon);
    uetype = unique(etype);
    for kk=1:length(uetype)
        imatch = find(strcmp(etype,uetype{kk}));
        stit = sprintf('%7i / %7i = %.3f are event type %s',length(imatch),n,length(imatch)/n,uetype{kk});
        disp(stit);
        figure; ax0 = [190 230 56 72];
        alaska_basemap(ax0,0,0,0)
        plot(wrapTo360(lon(imatch)),lat(imatch),'.');
        title(stit); axis(ax0);
    end
    % subset of quarry blasts
    iQ = find(strcmp(etype,'Q')==1);
    figure; plot(lon(iQ),lat(iQ),'.');
    
end

%==========================================================================
