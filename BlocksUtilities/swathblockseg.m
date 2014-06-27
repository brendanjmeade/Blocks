function [sx, sy] = swathblockseg(seg, b, d)
% SWATHSEG  Produces a polygon bounding a series of segmentseg.
%   [SX, SY] = SWATHBLOCKSEG(SEG, BOC, D) produces a polygon bounding the segments
%   in structure SEG, located +/- D km from the segmentseg.  At each segment
%   endpoint, a line is projected D km at +/- 90 degrees to the mean strike
%   of the 2 adjoining segmentseg.  No error checking is done to ensure non-
%   overlapping normal pointseg.  SEG is an ordered subset of segments 
%   comprising a single block, assuming that the first and last endpoints are
%   connected.  BOC is an n-by-2 array containing the ordered [lon. lat.] 
%   coordinates from fields "orderLon" and "orderLat" of the block structure
%   returned from BLOCKLABEL.
%

% Augment the segment array so that all segments have a "previous" segment
nseg = length(seg.lon1); % number of segments
slast = structsubset(seg, nseg);
seg = structmath(slast, seg, 'vertcat');

nseg = nseg + 1; % re-count number of segments

% Determine azimuths - using Cartesian instead of geographic
az = 90 - rad2deg(atan2((seg.lat2 - seg.lat1), (seg.lon2 - seg.lon1)));
%az = azimuth(seg.lat1(:), seg.lon1(:), seg.lat2(:), seg.lon2(:));

% Ordered azimuths
[bo, ba] = deal(b(:, 1), b(:, 2));
bo = [bo; bo(1)]; ba = [ba; ba(1)];
caz = wrapTo360(90 - rad2deg(atan2((ba(2:end) - ba(1:end-1)), (bo(2:end) - bo(1:end-1)))));
%caz = azimuth(ba(1:end-1), bo(1:end-1), ba(2:end), bo(2:end));
% Check to see if we need to flip this
if ~ismember(bo(3), [seg.lon1(1:3); seg.lon2(1:3)])
   caz = flipud(caz);
end
dew = [90 - caz, caz - 270];
[junk, mew] = min(abs(dew), [], 2);
mew = sub2ind(size(dew), (1:size(dew, 1))', mew);
mew = dew(mew);
mew = [mew(1); mew];

% Perpendicular azimuths
gt90 = az > 90;
paze = az + 90;
paze(gt90) = az(gt90) - 90;
pazw = az - 90;
pazw(gt90) = az(gt90) + 90;

% Allocate space
sx = nan(2*(nseg - 1) + 2, 1);
sy = sx;

% Convert d from km to degrees for reckoning
d = rad2deg(d./6371);
de = 10*d;


for i = 2:nseg
   idx1 = i-1;
   idx2 = length(sx) - (i - 1) + 1;
   % Reckon off-endpoint points

   % Start point
   [ey1es, ex1es] = reckon(seg.lat1(i), seg.lon1(i), d, paze(i));
   [ey1ws, ex1ws] = reckon(seg.lat1(i), seg.lon1(i), d, pazw(i));
   
   % End point
   [ey1ee, ex1ee] = reckon(seg.lat2(i), seg.lon2(i), d, paze(i));
   [ey1we, ex1we] = reckon(seg.lat2(i), seg.lon2(i), d, pazw(i));
   
   % Start point
   [ey2es, ex2es] = reckon(seg.lat1(i-1), seg.lon1(i-1), d, paze(i-1));
   [ey2ws, ex2ws] = reckon(seg.lat1(i-1), seg.lon1(i-1), d, pazw(i-1));
   
   % End point
   [ey2ee, ex2ee] = reckon(seg.lat2(i-1), seg.lon2(i-1), d, paze(i-1));
   [ey2we, ex2we] = reckon(seg.lat2(i-1), seg.lon2(i-1), d, pazw(i-1));

   % Find intersections
   [xiee, yiee] = pbisectfull([ex1es ey1es], [ex1ee ey1ee], [ex2es ey2es], [ex2ee ey2ee]);
   [xiww, yiww] = pbisectfull([ex1ws ey1ws], [ex1we ey1we], [ex2ws ey2ws], [ex2we ey2we]);

   [xiew, yiew] = pbisectfull([ex1es ey1es], [ex1ee ey1ee], [ex2ws ey2ws], [ex2we ey2we]);
   [xiwe, yiwe] = pbisectfull([ex1ws ey1ws], [ex1we ey1we], [ex2es ey2es], [ex2ee ey2ee]);
   
%   plot(xiee, yiee, '.r')
%   plot(xiww, yiww, '.m')
%
%   plot(xiew, yiew, '.g')
%   plot(xiwe, yiwe, '.c')   
   
   % Choose intersection points based on azimuths
   
   % Default case is that ee and ww are chosen; we'll do some resorting later
   [sx(idx1), sy(idx1)] = deal(xiee, yiee);
   [sx(idx2), sy(idx2)] = deal(xiww, yiww);

   % Change segments that flip sides of E-W to be other points
   if sign(mew(i)) ~= sign(mew(i-1)) 
      [sx(idx1), sy(idx1)] = deal(xiew, yiew);
      [sx(idx2), sy(idx2)] = deal(xiwe, yiwe);
   end
end

[sx(nseg), sy(nseg)] = deal(sx(1), sy(1));
[sx(nseg+1), sy(nseg+1)] = deal(sx(end), sy(end));

% Wrap to input convention
sx = wrapTo360(sx);

% Check for segment crossing; flip if need be
sx2 = sx; sy2 = sy;
for i = 2:length(seg.lon1)
   [ix(i), iy] = pbisect([sx(i-1) sy(i-1)], [sx(i) sy(i)], [seg.lon1(i) seg.lat1(i)], [seg.lon2(i) seg.lat2(i)]);
   if ~isnan(ix(i))
      sx(i) = sx2(2*nseg - i + 1);
      sy(i) = sy2(2*nseg - i + 1);
      sx(2*nseg - i + 1) = sx2(i);
      sy(2*nseg - i + 1) = sy2(i);
   end
end
