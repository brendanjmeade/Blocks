function h = smoothblockcline(seg, b, wid, col)
% SMOOTHBLOCKCLINE  Plots a smooth colored "line" as a series of polygons
%    SMOOTHBLOCKCLINE(SEG, B, WID, COL) uses the segment structure SEG and the 
%    ordered block coordinates B ([orderLon orderLat as returned by BLOCKLABEL)
%    to plot colored polygons of width WID (km), colored by the magnitude of
%    vector COL around a block.  For example, to plot a block colored by slip 
%    rate, 10 km wide, call:
%   
%    >> smoothcline(seg, [blon, blat], 10, seg.ssRate);
%
%    where SEG is a subset of the entire Mod.segment structure.  
%

% Define the polygons
[sx, sy] = swathblockseg(seg, b, wid/2);

% Define the individual polygons' indices
nseg = length(seg.lon1);
idx = [1:nseg; 2:nseg+1; fliplr(nseg+2:2*nseg+1); fliplr(nseg+3:2*nseg+2)];

% Make the plot
figure
h = patch('vertices', [sx(:) sy(:)], 'faces', idx', 'facevertexcdata', col);
shading flat;