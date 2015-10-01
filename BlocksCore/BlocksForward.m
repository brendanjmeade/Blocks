function v = BlocksForward(x, y, outdir, notri)
%
% BLOCKSFORWARD calculates the predicted velocities at stations for 
% a given BLOCKS model run.
%    BLOCKSFORWARD(X, Y, OUTDIR) calculates the predicted velocity 
%    components at station coordinates (X, Y) given the results of
%    a BLOCKS run contained in OUTDIR.
%
%    V = BLOCKSFORWARD(...) outputs the predicted velocities to structure
%    V with fields Mod, Rot, Def, and Tri. Each field is a 3*nStations-by-1
%    array arranged as [E1, N1, U1, ... En, Nn, Un]'
%

% Read in the necessary files
cmname                                           = dir([outdir filesep '*.command']);
Command                                          = ReadCommand([outdir filesep cmname(1).name]);
Segment                                          = ReadSegmentTri([outdir filesep 'Mod.segment']);
Segment                                          = ProcessSegment(Segment, Command);
Block                                            = ReadBlock([outdir filesep 'Mod.block']);
[Patches.c, Patches.v, Patches.s]                = PatchData([outdir filesep 'Mod.patch']);
Patches.nc                                       = size(Patches.c, 1);
Patches.nEl                                      = size(Patches.v, 1);
if Patches.nc ~= 0
   Patches                                       = PatchCoords(Patches);
end
Command                                          = [];

if exist('notri', 'var')
   Segment.patchTog = 0*Segment.patchTog;
end

% Put the stations in a structure
[Station.lon, Station.lat]                       = deal(x(:), y(:));
[Station.x, Station.y, Station.z]                = sph2cart(DegToRad(Station.lon), DegToRad(Station.lat), 6371);
[Segment, Block, Station]                        = BlockLabel(Segment, Block, Station);


% Calculate all partials
Partials.elastic                                 = GetElasticPartials(Segment, Station);
Partials.slip                                    = GetSlipPartials(Segment, Block);
Partials.rotation                                = GetRotationPartials(Segment, Station, Command, Block);
if sum(Segment.patchTog) > 0 % if patches are involved at all
   [Partials.tri, tz, ts]                        = GetTriCombinedPartials(Patches, Station, [1 0]);
   ts = zeros(3*size(Patches.s, 1), 1);
   ts(1:3:end) = Patches.s(:, 1);
   ts(2:3:end) = Patches.s(:, 2);
   ts(3:3:end) = Patches.s(:, 3);
else
   Partials.tri                                  = zeros(3*numel(x), 1);
   ts                                            = 0;
end

% Construct solution vectors
[x, y, z] = EulerToOmega(Block.eulerLon(:), Block.eulerLat(:), Block.rotationRate*1e6);
omega = zeros(3*numel(x), 1);
omega(1:3:end) = x(:);
omega(2:3:end) = y(:);
omega(3:3:end) = z(:);

% Do the forward problems
v.Rot = Partials.rotation*omega;
v.Def = Partials.elastic*Partials.slip*omega;
v.Tri = Partials.tri*ts;
v.Mod = v.Rot - v.Def - v.Tri; % Total modeled velocities