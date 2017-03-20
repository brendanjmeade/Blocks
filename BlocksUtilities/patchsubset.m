function pn = patchsubset(p, idx, except)
% PATCHSUBSET  Returns a subset of elements contained in patch structure.
%   PN = PATCHSUBSET(P, IDX) creates a new structure, PN, that contains
%   the subset of elements of structure P defined by indices IDX.  IDX 
%   can be a logical array of length SUM(P.nEl) or an array of indices into
%   P.v.
%

% Parse optional exceptions
if exist('except', 'var')
   ign = cat(2, {'c', 'nc', 'nEl', 'up', 'uue'}, except(:)');
else 
   ign = {'c', 'nc', 'nEl', 'up', 'uue'};
end

% Take the subset
pn = structsubset(p, idx, ign);

% Determine how many elements of each individual geometry are retained 
% (i.e., define p.nEl)

% Convert to real indices
if islogical(idx)
   idx = find(idx);
end

% Now trim coordinates
newc = pn.c(unique(pn.v), :);
if isfield(p, 'lon1')
   [tf1, newv(:, 1)] = ismember([pn.lon1 pn.lat1 pn.z1], newc, 'rows');
   [tf2, newv(:, 2)] = ismember([pn.lon2 pn.lat2 pn.z2], newc, 'rows');
   [tf3, newv(:, 3)] = ismember([pn.lon3 pn.lat3 pn.z3], newc, 'rows');
elseif isfield(p, 'x1')
   [tf1, newv(:, 1)] = ismember([pn.x1 pn.y1 pn.z1], newc, 'rows');
   [tf2, newv(:, 2)] = ismember([pn.x2 pn.y2 pn.z2], newc, 'rows');
   [tf3, newv(:, 3)] = ismember([pn.x3 pn.y3 pn.z3], newc, 'rows');
end

% Determine how many indices are less than the cumulative number of 
% elements for each geometry
cnel = cumsum([0; p.nEl(:)]);
for i = 1:numel(p.nEl)
   pn.nEl(i) = sum(idx <= cnel(i+1) & idx > cnel(i));
end

% Update number of element array to remove any 0 values
% This means that entire geometry was removed
if isfield(pn, 'up')
   pn.up = pn.up(pn.nEl > 0);
end
pn.nEl = pn.nEl(pn.nEl > 0);

pn.v = newv;
pn.c = newc;
ends = cumsum(pn.nEl);
begs = [1; ends(1:end-1) + 1];
pn.nc = zeros(size(pn.nEl));
for i = 1:numel(pn.nEl)
   pn.nc(i) = length(unique(pn.v(begs(i):ends(i), :)));
end
