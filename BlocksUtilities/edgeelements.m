function [top, bot, s1, s2] = edgeelements(c, v, strthresh)
% EDGEELEMENTS  Finds elements lining the edge of a mesh.
%   [top, bot, s1, s2] = EDGEELEMENTS(c, v, thresh) finds the 
%   indices of elements lining the edges of the mesh defined by 
%   coordinate array c and vertex ordering array v, returning the 
%   indices to top, bot, s1, and s2. This algorithm does not rely on 
%   finding the tops and bottoms based on nodes lying at a particular 
%   depth, instead relying on finding abrupt changes in the trends of 
%   the mesh edges, set by optional third argument thresh (the default
%   is 55 degrees). This makes it more suitable for meshes that have 
%   irregular depth top traces (such as subduction trenches) and/or are
%   not generated based on a set of depth contours.
%
% 

% Define default strike change threshold
if ~exist('strthresh', 'var')
   strthresh = 55;
end

% Allocate space for edge matrices
[top, bot, s1, s2] = deal([]);

% Get ordered edge coordinates
elo = OrderedEdges(c, v);
% Use ordered coordinates to calculate edge strikes
str = sphereazimuth(c(elo(1, :), 1), c(elo(1, :), 2), c(elo(2, :), 1), c(elo(2, :), 2));
elo = [elo(1, :) elo(1, 1)]; 
% Calculate difference in strikes of adjacent edges
dstr = diff([str; str(1)]);
% Find sharp angles; these are mesh corners. These are indices into elo.
corn = find(abs(dstr) > strthresh & abs(dstr) < (360-strthresh));
corn = [corn; corn(1)];

% Use corner indices to separate edges. Make inherent assumption that the top and bottom are longer than the sides
dcorn = diff(corn);
% Address negative index differences by adding size of elo
% Adjusted corners, just for differencing
corna = corn;
% Find the negative difference and add size of elo to all subsequent corner indices
corna(find(dcorn < 0, 1)+1:end) = corna(find(dcorn < 0, 1)+1:end) + size(elo, 2);
dcorn = diff(corna);
edgeend = cumsum(dcorn); edgeend(end) = edgeend(end) - 1;
edgebeg = [1; edgeend(1:end-1)];
[~, sidx] = sort(dcorn); % Get the indices of the longer edges

% Loop through and find the elements along each edge
for j = 1:length(dcorn)
   % Indices of all coordinates along this edge
   eidx = 1+(corn(j):corn(j+1));
   if isempty(eidx) % This happens when the second index is less than the first index
      eidx = 1+[(corn(j)):(length(elo)-1), 0:corn(j+1)];
   end
   % Find elements that contain at least 2 of these coordinates
   eedge = sum(ismember(v, elo(eidx)), 2) >= 2;
 
   % Assign to distinct edges
   if ismember(j, sidx(1:2)) % Working on a side 
      if isempty(s1)
         s1 = eedge;
      else
         s2 = eedge;
      end
   else
      % Working on a top or bottom, so let's test the depth of the corner
      if c(elo(corn(j+1)), 3) < mean(c(:, 3))
         bot = eedge;
      else
         top = eedge;
      end
   end
end