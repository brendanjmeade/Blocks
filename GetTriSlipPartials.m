function [G, Patches] = GetTriSlipPartials(Patches, Block, Segment)
% Calculate slip partial derivatives for triangular elements
nPatches                                         = sum(Patches.nEl);
nBlocks                                          = numel(Block.eulerLon);
G                                                = zeros(3*nPatches, 3*nBlocks);
[Patches.patchNum, ...     
 Patches.eastLabel, ...     
 Patches.westLabel]                              = deal(zeros(nPatches, 1));
[Patches.xxc, Patches.yyc, Patches.zzc]          = sph2cart(deg2rad(Patches.lonc), deg2rad(Patches.latc), 6371 - abs(Patches.zc));
 
% Determine which blocks lie on either side of each patch

% Find segments that are replaced by patches 
psi                                              = intersect(find(Segment.patchFile), find(Segment.patchTog));
Seg                                              = structsubset(Segment, psi);

% Switch the zeroing array
tz                                               = Patches.tz;
tz(Patches.tz == 2)                              = 3;
tz(Patches.tz == 3)                              = 2;

% Strike correction: Elements that dip > 90 
corr                                             = ones(nPatches, 1);
corr(Patches.strike > 270 | Patches.strike < 90) = 0;

for iPatch = 1:nPatches
   % Find which mesh this TDE belongs to
   Patches.patchNum(iPatch)                      = find(1+[0; cumsum(Patches.nEl)] <= iPatch, 1, 'last');
   % Find all segments corresponding to that mesh
   nseg                                          = find(Seg.patchFile == Patches.patchNum(iPatch));
   useLons                                       = Seg.midLon(nseg);
   useLats                                       = Seg.midLat(nseg);
   % Find the segment that is closest to this element
   nsegInd                                       = dsearchn([useLons(:) useLats(:)], [Patches.lonc(iPatch) Patches.latc(iPatch)]);
   % Assign segment and its block labels to a new field of the structure
   Patches.nearSeg(iPatch)                       = nseg(nsegInd);
   Patches.eastLabel(iPatch)                     = Seg.eastLabel(Patches.nearSeg(iPatch));
   Patches.westLabel(iPatch)                     = Seg.westLabel(Patches.nearSeg(iPatch));
   
   % Projection from Cartesian to spherical coordinates at element centroids
   rowIdx                                        = (iPatch-1)*3+1;
   colIdxE                                       = (Patches.eastLabel(iPatch)-1)*3+1;
   colIdxW                                       = (Patches.westLabel(iPatch)-1)*3+1;
   R                                             = GetCrossPartials([Patches.xxc(iPatch) Patches.yyc(iPatch) Patches.zzc(iPatch)]);
   [vn_wx, ve_wx, vu_wx]                         = CartVecToSphVec(R(1,1), R(2,1), R(3,1), Patches.lonc(iPatch), Patches.latc(iPatch));
   [vn_wy, ve_wy, vu_wy]                         = CartVecToSphVec(R(1,2), R(2,2), R(3,2), Patches.lonc(iPatch), Patches.latc(iPatch));
   [vn_wz, ve_wz, vu_wz]                         = CartVecToSphVec(R(1,3), R(2,3), R(3,3), Patches.lonc(iPatch), Patches.latc(iPatch));

   % Use strike and dip to resolve slip onto plane
   R                                             = [ve_wx vn_wx vu_wx; ve_wy vn_wy vu_wy; ve_wz vn_wz vu_wz]';
   R                                             = ProjectDispPartialsMats(R, Patches.strike(iPatch), Patches.dip(iPatch));
%    strike                                        = Patches.strike(iPatch) - 180.*corr(iPatch);
%    rot                                           = [sind(strike), cosd(strike), 0; cosd(strike), -sind(strike), 0; cosd(strike), -sind(strike), 0];
%    R                                             = rot*R;
   R(2, :)                                       = R(2, :)./cosd(Patches.dip(iPatch)); 
   R(tz(iPatch), :)                              = [0 0 0];

   % Place into appropriate columns
   G(rowIdx:rowIdx+2,colIdxE:colIdxE+2)          = R;
   G(rowIdx:rowIdx+2,colIdxW:colIdxW+2)          = -R;
end
