function epdiff(dirnames, refidx)
% epdiff  Makes a plot showing Euler pole differences. 
%    epdiff(DIRNAMES) makes a plot showing the differences
%    in estimated Euler poles contained in the Mod.block files
%    within Blocks run directories DIRNAMES. DIRNAMES is a 
%    character array giving the full path of the results directories
%    to be compared. By default, the difference with respect to the
%    first directory in the list is used as the reference directory
%    to which other runs' Euler poles are compared.
%
%    epdiff(DIRNAMES, REFIDX) allows specification of the row index
%    of the directory to be used as the reference as REFIDX.
%

% Check whether a reference index has been specified
if ~exist('refidx', 'var')
   refidx = 1;
end

ndirs = size(dirnames, 1);
otheridx = setdiff(1:ndirs, refidx);

% Load the reference block structure
b(refidx) = ReadBlock([dirnames(refidx, :) filesep 'Mod.block']);
nblocks = size(b(refidx).name, 1);
[pn, pe, pd, prot, angdiff, rotdiff, loc] = deal(zeros(nblocks, ndirs));
bc = ReadBlockCoords(dirnames(refidx, :));

[~, loc(:, refidx)] = 1:nblocks;
[pn(:, refidx), pe(:, refidx), pd(:, refidx)] = sph2cart(deg2rad(b(refidx).eulerLon), deg2rad(b(refidx).eulerLat), 1);
prot(:, refidx) = b(refidx).rotationRate;

figure; hold on
cmap = jet(ndirs-1);
j = length(otheridx);
% Load in the other block structures 
for i = fliplr(otheridx)
   b(i) = ReadBlock([dirnames(i, :) filesep 'Mod.block']);
   % Match blocks to reference by name
   [~, loc(:, i)] = ismember(b(refidx).name, b(i).name, 'rows');
   % Convert pole to Cartesian coordinates
   [pn(:, i), pe(:, i), pd(:, i)] = sph2cart(deg2rad(b(i).eulerLon(loc(:, i))), deg2rad(b(i).eulerLat(loc(:, i))), 1);
   prot(:, i) = b(i).rotationRate(loc(:, i));
   % Angular difference in Euler pole
   angdiff(:, i) = acosd(dot([pn(:, i), pe(:, i), pd(:, i)], [pn(:, refidx), pe(:, refidx), pd(:, refidx)], 2));
   angdiff(angdiff(:, i) > 90, i) = 180 - angdiff(angdiff(:, i) > 90, i);
   rotdiff(:, i) = prot(:, i) - prot(:, refidx);
   scatter(1:nblocks, angdiff(:, i), 5+ceil(50*abs(rotdiff(:, i))), cmap(j, :), 'filled', 'markeredgecolor', 'k');
   j = j-1;
end

cb = colorbar;
caxis([0 ndirs]); colormap(jet);

prepfigprint;
ylabel(cb, 'Run number');
set(cb, 'ticks', 1:ndirs-1);
xlabel('Block number')
ylabel('Angular difference in Euler pole (degrees)')
axis tight




