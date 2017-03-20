function pn = removeisectmesh(p)
% REMOVEISECTMESH  Removes small intersecting meshes.
%    pn = REMOVEISECTMESH(p) removes intersecting meshes
%    from mesh structure p. When two meshes intersect, 
%    whichever has fewer elements is removed. The trimmed
%    structure is returned to pn. 
%

% Find element intersections
crossidx = meshcrossing(p.c, p.v);

% Find which mesh each crossed element belongs to
meshidx = reshape(idmeshes(p, crossidx(:)), size(p.v));

% Find which mesh all elements belong to
allmeshidx = idmeshes(p, (1:size(p.v, 1))');

% Compare mesh sizes, discarding the smallest
toss = [];
for i = 1:size(meshidx, 1)
   if sum(meshidx(i, :), 2) > 0
      if p.nEl(allmeshidx(i)) < min(p.nEl(meshidx(i, meshidx(i, :) > 0)))
         toss = [toss; allmeshidx(i)];
      end
   end
end

% Trim the patch structure
pn = patchsubset(p, ~ismember(allmeshidx, toss));