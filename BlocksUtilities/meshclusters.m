function clus = meshclusters(p, sel)
% meshclusters   Finds isolated clusters of selected elements.
%   meshclusters(P, SEL) finds clusters of elements within the 
%   patch structure P that are selected, as identified in vector
%   SEL, which contains indices of selected elements or is a 
%   logical array giving selected elements. The clusters' element
%   indices are returned to a cell. 
%
%   If all elements are selected in SEL, then the returned indices
%   are those of all elements. 
%

% Blank cluster index array
idx = zeros(0, 2); lidx = [];

% Find the bounding edges of the selected elements
el = boundedges(p.c, p.v(sel, :));

% Check to see if there are any "necked" clusters (meeting at a single vertex)
uel = unique(el(:));
nn = hist(el(:), uel)';
necks = [uel(find(nn > 2)), nn(find(nn > 2))];
El = el;
elvec = 1:size(El, 1);
% Start with necks
if ~isempty(necks)
   for i = 1:size(necks, 1)
      start = find(el(:, 1) == necks(i, 1));
      for j = 1:length(start)
         elo = el(start(j), :);
         if ismember(elo, El, 'rows')
				for k = 2:elvec(end)
					[next, col] = find(El == elo(k-1, 2)); % find all of the boundary lines containing the second entry of the current ordered boundary line
					n = find(sum(El(next, :), 2) ~= sum(elo(k-1, :), 2)); % choose that which is not the current boundary line
					next = next(n(1)); col = col(n(1));
					if col == 1
						elo(k, :) = El(next, :);
					else
						elo(k, :) = El(next, [2 1]);
					end
				end
				% Find the unique edge indices. All this should do is to eliminate duplicate edge entries, i.e., circling around the cluster again
				[uelo, ix] = unique(elo, 'rows', 'stable');
				nnn(i, :) = [i sum(ismember(uelo, idx, 'rows')) size(uelo, 1)];
				% If not all of those unique edges are already included in the cluster array, then add them
				if sum(ismember(uelo, idx, 'rows')) ~= size(uelo, 1) 
					idx = [idx; uelo];
					lidx = [lidx; size(uelo, 1)];
				end
				% After each loop, setdiff to remove the used edges.
				El = setdiff(El, idx, 'rows');
				elvec = 1:size(El, 1);
         end
      end   
   end          
end

% Now handle all other edges, using a while loop until we use up all available edges
i = 1;
while size(El, 1) > 3 % The 4 ignores single-element holes
   elo = El(1, :);
   for j = 2:elvec(end)
      [next, col] = find(El == elo(j-1, 2)); % find all of the boundary lines containing the second entry of the current ordered boundary line
      n = find(sum(El(next, :), 2) ~= sum(elo(j-1, :), 2)); % choose that which is not the current boundary line
      next = next(n(1)); col = col(n(1));
      if col == 1
         elo(j, :) = El(next, :);
      else
         elo(j, :) = El(next, [2 1]);
      end
   end
   % Find the unique edge indices. All this should do is to eliminate duplicate edge entries, i.e., circling around the cluster again
   [uelo, ix] = unique(elo, 'rows', 'stable');
   nnn(i, :) = [i sum(ismember(uelo, idx, 'rows')) size(uelo, 1)];
   % If not all of those unique edges are already included in the cluster array, then add them
   if sum(ismember(uelo, idx, 'rows')) ~= size(uelo, 1) 
      idx = [idx; uelo];
      lidx = [lidx; size(uelo, 1)];
   end
   % After each loop, setdiff to remove the used edges.
   El = setdiff(El, idx, 'rows');
   elvec = 1:size(El, 1);
   i = i+1;
end

% Now loop through all cluster edges and find elements lying within
last = cumsum(lidx);
first = [1; last(1:end-1)+1];
if islogical(sel)
   sel = find(sel);
end
j = 1;
for i = 1:length(lidx)
   clu = intersect(sel, find(inpolygon(p.lonc, p.latc, p.c(idx(first(i):last(i), :), 1), p.c(idx(first(i):last(i), :), 2))));
   if ~isempty(clu)
      clus{j} = clu;
      j = j+1;
   end
end

