function S = snapsegments(s, p, sel, dtol)
% snapsegments  Snaps segments to trace the updip edge of patches.
%   snapsegments(S, P, SEL, DTOL) moves segments in structure S so that
%   they coincide with the updip edge of a patch described in structure P.
%   Only the selected segments with indices SEL will be altered. Those 
%   segments are removed from the structure and replaced with the actual 
%   triangle edge coordinates from P. DTOL is an optional depth tolerance
%   variable and is used to allow some variation in choosing the updip edge
%   of the meshes. DTOL is a constant specifying the depth tolerance in km. 
%
%   The function operates by finding the patch whose centroid is closest to 
%   the mean of the selected segment coordinates. Then, the selected segments
%   are split or removed and replaced with new segments that overlap with the
%   updip trace of the patch. These segments have "blank" properties, i.e., 
%   the names, locking depths, dips, etc. are those of default segments, which
%   should not matter in a Blocks run since these segments will not contribute
%   to the elastic deformation field.
%
%   When selecting the segments that will be replaced by the patches, it is best
%   to select those that fully encompass the updip extent of the mesh. That is,
%   make sure that the selected segments' endpoints lie beyond (even slightly) the
%   extent of the mesh.
%

% Find center of selected segments
sc = [mean([s.lon1(sel); s.lon2(sel)]), mean([s.lat1(sel); s.lat2(sel)])];

% Find centroids of meshes
mn = length(p.nEl);
   % Start and end indices
vends = cumsum(p.nEl(:));
vbegs = [1; vends(1:end-1)+1];

mc = zeros(mn, 2);
for i = 1:mn
   mc(i, :) = mean(p.c(p.v(vbegs(i):vends(i), :), 1:2), 1);
end

% Find the mesh closest to the segments
dist = gcdist(sc(2), sc(1), mc(:, 2), mc(:, 1));
[~, midx] = min(dist);



% Find the updip elements of the selected mesh
%cends = cumsum(p.nc(:));
%cbegs = [1; cends(1:end-1)+1];
%elRange                                   = vbegs(midx):vends(midx);
%coRange                                   = cbegs(midx):cends(midx);
%% find the zero depth coordinates
if ~exist('dtol', 'var')
   dtol = 0;
end
%zeroTest                               = abs(p.c(coRange, 3)) <= dtol;
%zeroZ                                  = cbegs(midx) + find(zeroTest);
%% find those elements having two zero-depth coordinates
%updip                                  = vbegs(midx) + find(sum(ismember(p.v(elRange, :), zeroZ), 2) == 2);

% Find the segments that are at the ends of the selection. These will be split, with the new endpoint
% coincident with the corner of the mesh.

% All coordinates
lons										= [s.lon1(sel(:)); s.lon2(sel(:))];
lats = [s.lat1(sel(:)); s.lat2(sel(:))];
coords									= [lons(:) lats(:)];
% Unique coordinates
[~, ucInd1]										= unique(coords, 'rows', 'first');
[uc, ucInd2]										= unique(coords, 'rows', 'last');
% Dangling segments are those where the unique occurrence is in the same place
[endsegs, endcol] = ismember(coords, uc(ucInd1 == ucInd2, :), 'rows');
endsegs = find(endsegs);
endsegs(endsegs > length(sel)) = endsegs(endsegs > length(sel)) - length(sel);
endcol = endcol(endcol > 0); % Tells whether the hanging point is endpoint 1 or 2
lons = reshape(lons, length(sel), 2);
lats = reshape(lats, length(sel), 2);

% Find corner coordinates of mesh. 
% - Find ordered edges
% - Find those along updip edge
% - Use distance to associate the first and last with the dangling segments
elo = OrderedEdges(p.c, p.v(vbegs(midx):vends(midx), :));
elo = elo(1, [2:end, 1]); % Ordered edge nodes
updip = elo(find(abs(p.c(elo, 3)) <= dtol)); % Updip nodes

% Find the corners
d1 = gcdist(lats(endsegs(1), endcol(1)), lons(endsegs(1), endcol(1)), p.c(updip, 2), p.c(updip, 1));
[~, c1] = min(d1);
d2 = gcdist(lats(endsegs(2), endcol(2)), lons(endsegs(2), endcol(2)), p.c(updip, 2), p.c(updip, 1));
[~, c2] = min(d2);

% Split the hanging segments
% New endpoints are mesh corners
newseg.lon1 = zeros(length(updip)+1, 1);
newseg.lon2 = zeros(length(updip)+1, 1);
newseg.lat1 = zeros(length(updip)+1, 1);
newseg.lat2 = zeros(length(updip)+1, 1);

newseg.lat2(1) = p.c(updip(c1), 2);
newseg.lon2(1) = p.c(updip(c1), 1);
% Other endpoint is the hanging endpoint
newseg.lat1(1) = lats(endsegs(1), endcol(1));
newseg.lon1(1) = lons(endsegs(1), endcol(1));
newseg.lat2(2) = p.c(updip(c2), 2);
newseg.lon2(2) = p.c(updip(c2), 1);
newseg.lat1(2) = lats(endsegs(2), endcol(2));
newseg.lon1(2) = lons(endsegs(2), endcol(2));

% Remaining new segments are the mesh edges
newseg.lon1(3:end) = p.c(updip(1:end-1), 1);
newseg.lat1(3:end) = p.c(updip(1:end-1), 2);
newseg.lon2(3:end) = p.c(updip(2:end), 1);
newseg.lat2(3:end) = p.c(updip(2:end), 2);
newseg.name = sprintf('Patch%g_%g', [midx*ones(1, length(updip)+1); 1:length(updip)+1]);
% Stitch new and old segments together
old = structsubset(s, setdiff(1:length(s.lon1), sel));
S = AddGenericSegment(old, newseg.name, newseg.lon1, newseg.lat1, newseg.lon2, newseg.lat2);






   