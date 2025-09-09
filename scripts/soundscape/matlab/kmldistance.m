% Load the kml files as tables
if bproposal
    addpath(ddir);  % for kml2struct.m
    %stag = 'dead_mouse'; N = 94;
    stag = 'parks_hwy'; N = 306;
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

% Load a KML/KMZ file with kmz2struct function
filename1 = 'Alaska_Railroad.kml'; % The full or relative path to KML/KMZ file you wish to load
S1 = kmz2struct(filename1); % The data loaded from the KML/KMZ file in MATLAB structure format

% If a MATLAB table is preferred, use the MATLAB struct2table function
T1 = struct2table(S1);

% Display the table
%disp(T1)

filename2 = 'ZE_NODAL.kml'; % The full or relative path to KML/KMZ file you wish to load
S2 = kml2struct(filename2); % The data loaded from the KML/KMZ file in MATLAB structure format

% If a MATLAB table is preferred, use the MATLAB struct2table function
T2 = struct2table(S2);

% Display the table
%disp(T1)

% Extract the coordinates as matrices
P1 = [T1.Latitude, T1.Longitude];
P2 = [T2.Latitude, T2.Longitude];

% Initialize a matrix to store the distances
D = zeros(size(P1,1), size(P2,1));

% Loop over each pair of points and compute the distance
for i = 1:size(P1,1)
    for j = 1:size(P2,1)
        % Use the haversine formula to calculate the distance in km
        D(i,j) = haversine(P1(i,:), P2(j,:));
    end
end

% Find the minimum distance and its indices
[min_dist, ind] = min(D(:));
[row, col] = ind2sub(size(D), ind);

% Display the result
fprintf('The shortest distance is %.3f km between point %d in file1 and point %d in file2.\n', min_dist, row, col);
