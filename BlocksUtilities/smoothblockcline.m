function [Plon, Plat] = smoothblockcline(direc, wid, col, bidx)
% SMOOTHBLOCKCLINE  Plots smooth colored "lines" as series of polygons.
%    SMOOTHBLOCKCLINE(DIREC, WID, COL) uses information in the results directory
%    DIREC to plot colored polygons of width WID (km), colored by 
%    COL around a block. COL can be a valid field from a segment structure,
%    or a nSegments-by-1 vector. For example, to plot segments colored by strike 
%    slip rate, 10 km wide, call:
%   
%    >> smoothblockcline('0000000001', 10, 'ssRate');
%
%    SMOOTHBLOCKCLINE(DIREC, WID, COL, IDX) plots only the blocks specified
%    by vector IDX. By default, the exterior block is not plotted, as its 
%    coordinates are not included in Block.coords. 
%
%    [PLON, PLAT] = SMOOTHBLOCKCLINE(...) returns longitude and latitude 
%    coordinates to the 4-by-nSegments arrays PLON and PLAT, which can 
%    be plotted using:
%
%    patch(PLON, PLAT, color_vector)
%
%    where color_vector is any 1-by-nSegments vector. 
% 

% Read in necessary files
seg = ReadSegmentTri([direc '/Mod.segment']);
lab = ReadSegmentTri([direc '/Label.segment']); lab.eastBlock = lab.ssRate; lab.westBlock = lab.ssRateSig;
b   = ReadBlockCoords([direc '/Block.coords']);

% Define subset of blocks
if ~exist('bidx', 'var')
   bidx = 1:length(b);
end

% Calculate segment midpoint coordinates   
mlon = (seg.lon1 + seg.lon2)/2;   
mlat = (seg.lat1 + seg.lat2)/2;

% Big arrays to hold all polygon coordinates 
[Plon, Plat] = deal(zeros(4, length(seg.lon1)));
Daz = zeros(1, length(seg.lon1));
idaz = 0;

% Define the polygons for each block
for i = bidx
   [plon, plat] = swathblockseg(b{i}, wid);
   % Match the polygon coordinates to segments

   % Calculate ordered block midpoints
   bmid = (b{i} + b{i}([2:end, 1], :))/2;
   % Find coordinates within ordered block coordinates
   [~, loc] = ismember(bmid, [mlon, mlat], 'rows');
   % Reorder columns of polygon coordinates
   Plon(:, loc) = plon;
   Plat(:, loc) = plat;
end

% Correct any negative longitudes
Plon(Plon < 0) = Plon(Plon < 0) + 360;

% Make the plot

% Determine color variable
if ischar(col)
   if isfield(seg, col)
      col = getfield(seg, col);
   end
end
figure
h = patch(Plon, Plat, col(:)');
shading flat;