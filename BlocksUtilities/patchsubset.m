function pn = patchsubset(p, idx)
% PATCHSUBSET  Returns a subset of elements contained in patch structure.
%   PN = PATCHSUBSET(P, IDX) creates a new structure, PN, that contains
%   the subset of elements of structure P defined by indices IDX.  IDX 
%   can be a logical array of length SUM(P.nEl) or an array of indices into
%   P.v.
%

% Take the subset
pn = structsubset(p, idx, {'c', 'nc', 'nEl', 'up'});

% Determine how many elements of each individual geometry are retained 
% (i.e., define p.nEl)

% Convert to real indices
if islogical(idx)
   idx = find(idx);
end

% Determine how many indices are less than the cumulative number of 
% elements for each geometry
cnel = cumsum([0; p.nEl(:)]);
for i = 1:numel(p.nEl)
   pn.nEl(i) = sum(idx <= cnel(i+1) & idx > cnel(i));
end

pn.nc = sum(pn.nc);
pn.nEl = sum(pn.nEl);