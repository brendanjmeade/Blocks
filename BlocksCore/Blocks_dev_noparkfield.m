function Blocks(commandFile, trunc)

fprintf('Parsing input data...')
runName                                          = GetRunName; % Create new directory for output
Command                                          = ReadCommand(commandFile); % Read command file
Station                                          = ReadStation(Command.staFileName); % Read station file
if strmatch(Command.unitSigmas, 'yes')
   [Station.eastSig, Station.northSig]           = ones(numel(Station.lon), 1);
end
Segment                                          = ReadSegmentTri(Command.segFileName); % Read segment file
Segment.bDep                                     = zeros(size(Segment.bDep));
Block                                            = ReadBlock(Command.blockFileName); % Read block file
Segment                                          = OrderEndpoints(Segment); % Reorder segment endpoints in a consistent fashion
Patches                                          = struct('c', [], 'v', []);
Segment.patchTog(Segment.patchFile == 1)         = 0;
Segment.patchFile(Segment.patchFile == 1)        = 0;
Segment.patchFile(Segment.patchFile == 2)        = 1;
if ~isempty(Command.patchFileNames)
   Patches                                       = ReadPatches(Command.patchFileNames); % Read triangulated patch files
   Patches                                       = PatchEndAdjust(Patches, Segment); % Adjust patch end coordinates to agree with segment end points
   [Patches, Command]                            = TogglePatches(Patches, Segment, Command); % Adjust structure so that only the toggled-on patch data are retained
   Patches                                       = PatchCoords(Patches); % Create patch coordinate arrays
   if numel(Command.triSmooth) == 1
      Command.triSmooth                          = repmat(Command.triSmooth, 1, numel(Patches.nEl));
   elseif numel(Command.triSmooth) ~= numel(Patches.nEl)
      error('BLOCKS:SmoothNEqPatches', 'Smoothing magnitude must be a constant or array equal in size to the number of patches.');
   end   
end
[Segment.x1 Segment.y1 Segment.z1]               = sph2cart(DegToRad(Segment.lon1(:)), DegToRad(Segment.lat1(:)), 6371);
[Segment.x2 Segment.y2 Segment.z2]               = sph2cart(DegToRad(Segment.lon2(:)), DegToRad(Segment.lat2(:)), 6371);
[Segment.midLon Segment.midLat]                  = deal((Segment.lon1+Segment.lon2)/2, (Segment.lat1+Segment.lat2)/2);
[Segment.midX Segment.midY Segment.midZ]         = sph2cart(DegToRad(Segment.midLon), DegToRad(Segment.midLat), 6371);
Segment.lDep                                     = LockingDepthManager(Segment.lDep, Segment.lDepSig, Segment.lDepTog, Segment.name, Command.ldTog2, Command.ldTog3, Command.ldTog4, Command.ldTog5, Command.ldOvTog, Command.ldOvValue);
Segment.lDep                                     = PatchLDtoggle(Segment.lDep, Segment.patchFile, Segment.patchTog); % Set locking depth to zero on segments that are associated with patches
Segment                                          = SegCentroid(Segment);
Station                                          = SelectStation(Station);
[Station.x Station.y Station.z]                  = sph2cart(DegToRad(Station.lon), DegToRad(Station.lat), 6371);
fprintf('done.\n')

% Assign block labels and put sites on the correct blocks
fprintf('Labeling blocks...')
[Segment, Block, Station]                        = BlockLabel(Segment, Block, Station);
fprintf('...done.\n')

fprintf('Applying a priori block motion constraints...')
Partials.blockCon                                = GetBlockConstraintPartials(Block);
blockConIdx                                      = find(Block.aprioriTog);
blockCons                                        = zeros(size(Partials.blockCon, 1), 1);
blockConSigs                                     = zeros(size(Partials.blockCon, 1), 1);
if sum(blockConIdx) > 0
   [apLons, apLats, apRates]                        = deal(Block.eulerLon(blockConIdx), Block.eulerLat(blockConIdx), Block.rotationRate(blockConIdx).*1e6);
   [apbx, apby, apbz]                               = EulerToOmega(apLons, apLats, apRates);
   blockCons(1:3:end)                               = apbx;
   blockCons(2:3:end)                               = apby;
   blockCons(3:3:end)                               = apbz;
   [apbsx, apbsy, apbsz]                            = epoles_cov_to_omega_cov(apbx, apby, apbz, diag([deg_to_rad(Block.eulerLatSig(blockConIdx)), deg_to_rad(Block.eulerLonSig(blockConIdx)), deg_to_rad(Block.rotationRateSig(blockConIdx))]));
   blockConSigs(1:3:end)                            = apbsx;
   blockConSigs(2:3:end)                            = apbsy;
   blockConSigs(3:3:end)                            = apbsz;
end   
fprintf('done.\n')

fprintf('Applying a priori slip constraints...')
% Build a priori slip rate constraints
for i = 1:numel(Segment.lon1)
   if Segment.ssRateTog(i)==1
      fprintf(1, 'Strike-slip constraint   : rate=%6.2f, sigma=%6.2f %s\n', Segment.ssRate(i), Segment.ssRateSig(i), Segment.name(i,:));
   end
   if Segment.dsRateTog(i)==1
      fprintf(1, 'Dip-slip constraint      : rate=%6.2f, sigma=%6.2f %s\n', Segment.dsRate(i), Segment.dsRateSig(i), Segment.name(i,:));
   end
   if Segment.tsRateTog(i)==1
      fprintf(1, 'Tensile-slip constraint  : rate=%6.2f, sigma=%6.2f %s\n', Segment.tsRate(i), Segment.tsRateSig(i), Segment.name(i,:));
   end
end
Partials.slipCon                                 = GetSlipPartials(Segment, Block);
slipToggle                                       = zeros(1,size(Partials.slipCon, 1));
slipToggle(1:3:end)                              = Segment.ssRateTog(:);
slipToggle(2:3:end)                              = Segment.dsRateTog(:);
slipToggle(3:3:end)                              = Segment.tsRateTog(:);
slipToggle                                       = slipToggle(:);
slipConIdx                                       = find(slipToggle==1);
slipRates                                        = zeros(1,size(Partials.slipCon, 1));
slipRates(1:3:end)                               = Segment.ssRate(:);
slipRates(2:3:end)                               = Segment.dsRate(:);
slipRates(3:3:end)                               = Segment.tsRate(:);
slipRates                                        = slipRates(:);
slipSigs                                         = zeros(1,size(Partials.slipCon, 1));
slipSigs(1:3:end)                                = Segment.ssRateSig(:);
slipSigs(2:3:end)                                = Segment.dsRateSig(:);
slipSigs(3:3:end)                                = Segment.tsRateSig(:);
slipSigs                                         = slipSigs(:);
Partials.slipCon                                 = Partials.slipCon(slipConIdx, :);
fprintf('done.\n')

fprintf('Calculating design matrix components...')
% Check to see if exisiting elastic kernels can be used
[Partials.elastic, Partials.tri,...
 Partials.trizeros, Partials.tristrikes]         = CheckExistingKernels(Command, Segment, Station, Patches);

% Get Partial derivatives for elastic calculation
if isempty(Partials.elastic)
   fprintf('\n  Calculating elastic partials...')
   Partials.elastic                              = GetElasticPartials(Segment, Station);
end
szelastic                                        = size(Partials.elastic);

% Check whether or not we are saving the partials
if strcmp(Command.saveKernels, 'yes') == 1;
   save('tempkernels.mat', '-struct', 'Partials', '-mat');
end

fprintf('\n  Calculating slip partials...')
Partials.slip                                    = GetSlipPartials(Segment, Block);
szslip                                           = size(Partials.slip);
fprintf('\n  Calculating rotation partials...')
Partials.rotation                                = GetRotationPartials(Segment, Station, Command, Block);
szrot                                            = size(Partials.rotation);
rowkeep                                          = setdiff(1:szrot(1), [3:3:szrot(1)]);

if sum(Segment.patchTog) > 0 & ~isempty(Patches.c) % if patches are involved at all
   if isempty(Partials.tri)
      fprintf('\n  Calculating triangular partials...')
      [Partials.tri,...
       Partials.trizeros,...
       Partials.tristrikes]                      = GetTriPartials(Patches, Station);
      % Check whether or not we are saving the partials
      if strcmp(Command.saveKernels, 'yes') == 1;
         save('tempkernels.mat', '-struct', 'Partials', 'tri', 'trizeros', 'tristrikes', '-mat', '-append');
      end
   end

   % Calculate triangular slip partials
%   Partials.trislip                              = GetTriSlipPartials(Patches, Block, Segment, Partials.tristrikes);

   % determine constraints placed on triangular slips by segment slips
   [Wconss, Wconst]                              = deal(zeros(0, szrot(2)), zeros(0, size(Partials.tri, 2)));
   if Command.triKinCons == 1
      [Wconss, Wconst]                           = PatchSlipConstraint(Partials.slip, Segment, Patches);
   end   
   % set slip partials to zero on segments that are replaced by patches
   %Partials.slip                                = ZeroSlipPartials(Segment.patchFile, Segment.patchTog, Partials.slip);
   szslip                                        = size(Partials.slip);
   
   % adjust triangular partials
   triS                                          = [1:length(Partials.trizeros)]';
   triD                                          = find(Partials.trizeros(:) == 2);
   triT                                          = find(Partials.trizeros(:) == 3);
   colkeep                                       = setdiff(1:size(Partials.tri, 2), [3*triD-0; 3*triT-1]);
   Partials.tri                                  = Partials.tri(:, colkeep); % eliminate the partials that equal zero
   % set up and downdip edges to zero if requested
   if sum(Command.triEdge) ~= 0
      Ztri                                       = ZeroTriEdges(Patches, Command);
   else
      Ztri                                       = zeros(0, size(Partials.tri, 2));
   end
   % adjust the matrix so that constant rake is applied
   if ~isempty(Command.triRake)
      Partials.tri                               = RakeTriPartials(Partials.tri, Partials.tristrikes, Command.triRake, triD, triT);
   end
   
   fprintf('\n  Making the smoothing matrix...')
   % Make the triangular smoothing matrix
   share                                         = SideShare(Patches.v);
   dists                                         = TriDistCalc(share, Patches.xc, Patches.yc, Patches.zc); % distance is calculated in km
%   Wtri                                          = MakeTriSmooth(share, dists);
   [Wtri, Wseg, beta]                            = SmoothEdges(share, Patches, Segment, Partials.slip, Command);                 
   Wtri                                          = Wtri(:, colkeep); % eliminate the entries corresponding to empty partials
   Wtri                                          = Wtri(colkeep, :); % eliminate the entries corresponding to empty partials
   Wseg(3:3:end, :)                              = [];
   beta(3:3:end)                                 = [];
   Wconst                                        = Wconst(:, colkeep);
   Wcons                                         = [Wconss Wconst];
   Wcons(3:3:end, :)                             = [];
%   beta                                          = zeros(2*length(Patches.v), 1);
%   s                                             = [0 cumsum(Patches.nEl)];
%   for i = 1:numel(Patches.nEl)
%      dist                                       = dists(s(i)+1:s(i+1), :);
%      distScale                                  = mean(dist(find(dist)));
%      beta(2*(s(i)+1)-1:2*s(i+1))                = distScale^2*Command.triSmooth(i);
%   end
%   [Wseg, Wtri]                                  = deal(repmat(beta', size(Wseg, 1), 1).*Wseg, repmat(beta', size(Wtri, 1), 1).*Wtri);
else
   [Wseg, Wtri, Wcons, beta, Ztri]               = deal([]);
end
sztri                                            = size(Partials.tri);

% Calculate strain partials based on the method specified in the command file
fprintf('\n  Calculating strain partials...')
[Partials.strain, strainBlockIdx]                = deal([]);
[Model.lonStrain, Model.latStrain]               = deal(zeros(size(Block.interiorLon)));
switch Command.strainMethod
   case 1 % use the block centroid
      [Partials.strain, strainBlockIdx,...
       Model.lonStrain, Model.latStrain]         = GetStrainCentroidPartials(Block, Station, Segment);
   case 2 % solve for the reference coordinates
      [Partials.strain, strainBlockIdx]          = GetStrain56Partials(Block, Station, Segment);
   case 3 % solve for the reference latitude only
      [Partials.strain, strainBlockIdx]          = GetStrain34Partials(Block, Station, Segment);
   case 4 % use the simulated annealing approach to solve for reference coordinates inside the block boundaries
      [Partials.strain, strainBlockIdx]          = GetStrainSearchPartials(Block, Station, Segment);
end      

fprintf('\n  Assembling the Jacobian...')
% Assemble Jacobian

% First determine submatrix sizes
rsta                                             = 2*numel(Station.lon); % station rows
rbcons                                           = size(Partials.blockCon, 1); % Block constraint rows
rscons                                           = size(Partials.slipCon, 1); % Slip constraint rows
rtriw                                            = sztri(2); % Triangular smoothing rows
rtriz                                            = size(Ztri, 1); % Triangular edge constraint rows

cblock                                           = szslip(2); % Block/slip columns
ctri                                             = rtriw; % Triangle columns
cstrain                                          = size(Partials.strain, 2); % Strain columns

% Determine indices
rnum                                             = [rsta rbcons rscons rtriw rtriz];
cnum                                             = [cblock ctri cstrain];
ridx                                             = cumsum([0 rnum]);
cidx                                             = cumsum([0 cnum]);

for i = 1:5
   for j = 1:3
      rows{i, j}                                 = ridx(i) + (1:rnum(i));
      cols{i, j}                                 = cidx(j) + (1:cnum(j));
   end
end

% Allocate space
R                                                = zeros(ridx(end), cidx(end));

Rb                                               = Partials.rotation - Partials.elastic * Partials.slip;
Rb                                               = Rb(rowkeep, :);
Rt                                               = -Partials.tri(rowkeep, :);

R(rows{1, 1}, cols{1, 1})                        = Rb;
R(rows{1, 2}, cols{1, 2})                        = Rt;
R(rows{1, 3}, cols{1, 3})                        = Partials.strain;
R(rows{2, 1}, cols{2, 1})                        = Partials.blockCon;
R(rows{3, 1}, cols{3, 1})                        = Partials.slipCon;
R(rows{4, 2}, cols{4, 2})                        = Wtri;
R(rows{5, 2}, cols{5, 2})                        = Ztri;

%R                                                = Partials.rotation - Partials.elastic * Partials.slip;
%R                                                = [R -Partials.tri]; % Insert tri. partials ** Negative because they're just like the elastic partials
%R                                                = R(rowkeep, :); % Eliminate the vertical velocity components
%R                                                = [R ; [Partials.blockCon sparse(size(Partials.blockCon, 1), sztri(2))]]; % Add block motion constraints
%R                                                = [R ; [Partials.slipCon sparse(size(Partials.slipCon, 1), sztri(2))]]; % Add slip rate constraints
%R                                                = [R ; [Wseg, Wtri]]; % Add triangular smoothing matrix
%R                                                = [R ; Wcons]; % Add kinematic constraints
%R                                                = [R ; [sparse(size(Ztri, 1), szrot(2)) Ztri]]; % Add edge slip constraints
%strainPadding                                    = zeros(size(R, 1)-size(Partials.strain, 1), size(Partials.strain, 2));
%R                                                = [R [Partials.strain;strainPadding]];
fprintf('done.\n')

fprintf('Building the data vector...')
% Build data vector and weighting matrix
d                                                = zeros(3*numel(Station.lon), 1);
d(1:3:end)                                       = Station.eastVel;
d(2:3:end)                                       = Station.northVel;
d(3:3:end)                                       = Station.upVel;
d                                                = d(rowkeep); % Eliminate the vertical velocity components
d                                                = [d ; blockCons]; % Add block motion constraints
d                                                = [d ; slipRates(slipConIdx)]; % Add slip rate constraints
d                                                = [d ; sparse(size(Wtri, 1), 1)]; % Add zeros for tri. mesh Laplacian
d                                                = [d ; sparse(size(Wcons, 1), 1)]; % Add zeros for rect. constraints
d                                                = [d ; sparse(size(Ztri, 1), 1)]; % Add zeros for edge slip constraints

w                                                = zeros(3*numel(Station.lon), 1);
w(1:3:end)                                       = 1./Station.eastSig.^2;
w(2:3:end)                                       = 1./Station.northSig.^2;
w(3:3:end)                                       = 1./Station.upSig.^2;
w                                                = w(:);
w                                                = w(rowkeep);
w                                                = [w ; Command.blockConWgt.*1./blockConSigs.^2]; % Add block motion constraints
w                                                = [w ; Command.slipConWgt.*1./(slipSigs(slipConIdx)).^2]; % Add slip rate constraints
w                                                = [w ; beta]; % add the triangular smoothing vector
w                                                = [w ; ones(size(Wcons, 1), 1)]; % add the triangular kinematic consistency vector
w                                                = [w ; 1e5*ones(size(Ztri, 1), 1)]; % add the zero edge vector
W                                                = spdiags(w, 0, numel(w), numel(w)); % assemble into a matrix
fprintf('done.\n')

% Estimate rotation vectors, Euler poles and slip rates
fprintf('Doing the inversion...\n')
%Command.triCons = [0 150];
Command.triCons = [];

% Check memory usage - see if we can make sparse matrices full
s1 = whos;
new = sum([size(R, 2)^2 ... % Model covariance
           size(Partials.slip, 2) size(Partials.tri, 2) size(Partials.strain, 2) ... % omegaEsts
           size(Partials.slip, 2)^2 size(Partials.tri, 2)^2 size(Partials.strain, 2)^2 ... % covariances
           5*size(Partials.slip, 2) ... % Euler poles
           2*size(Partials.slip, 1) ... % slip rates
           7*size(Partials.elastic, 1) ... % velocities
           3*size(Partials.tri, 2) ... % triangle slips
           10/3*size(Partials.strain, 2) ... % strain tensor parameters
           2*size(Partials.rotation, 2)]); % block strain rates
s1 = checkmemory(s1, new);
s2 = checkmemory(whos('R', 'W', 'd')); % find memory used by just the big sparse matrices
s2 = (s1 - s2 + 8*sum([numel(R) numel(W) numel(d)]))/1024^3; % this is the memory used with full matrices
fprintf('  Total memory usage is %3.2f Gb - using ', s2);
if s2 < 16
   fprintf('full matrices...\n')
   [R, W, d]                                     = deal(full(R), full(W), full(d));
else
   fprintf('sparse matrices...\n')
end


if isempty(Command.triCons)
   if ~exist('trunc', 'var')
      Model.covariance                           = (R'*W*R)\eye(size(R, 2));
      omegaEst                                   = Model.covariance*R'*W*d;
      keyboard
   else   
      [U, L, V]                                  = svd(sqrt(W)*R);
      Up                                         = U(:, 1:trunc);
      Lp                                         = L(1:trunc, 1:trunc); 
      Vp                                         = V(:, 1:trunc);
      iLp                                        = Lp\eye(length(Lp));
      Model.covariance                           = Vp*iLp.^2*Vp';
      omegaEst                                   = Vp*iLp*Up'*(sqrt(W)*d);
      keyboard
   end
else
   [minconstr, maxconstr]                        = deal(-inf(size(R, 2), 1), inf(size(R, 2), 1));
   minconstr(szslip(2)+1:szslip(2)+sztri(2))     = Command.triCons(1);
   maxconstr(szslip(2)+1:szslip(2)+sztri(2))     = Command.triCons(2);
   [omegaEst, resnorm, res, exf, output, lambda] = lsqlin(R, d, [], [], [], [], minconstr, maxconstr);
end
fprintf('\nInversion is done.\n')

fprintf('Calculating model results...')
% extract different parts of the state vector
omegaEstRot                                      = omegaEst(1:szslip(2)); % rotation parameters
omegaEstTriSlip                                  = omegaEst(szslip(2) + (1:sztri(2))); % triangular slips
omegaEstStrain                                   = omegaEst(szslip(2) + sztri(2) + 1:end); % strain parameters

% extract different parts of the model covariance matrix
Model.covarianceRot                              = Model.covariance(1:szslip(2), 1:szslip(2));
Model.covarianceTriSlip                          = Model.covariance(szslip(2) + (1:sztri(2)), szslip(2) + (1:sztri(2)));
Model.covarianceStrain                           = Model.covariance(szslip(2) + sztri(2) + 1:end, szslip(2) + sztri(2) + 1:end);

% Calculate Euler poles and uncertainties
[Model.rateEuler, Model.lonEuler, Model.latEuler]= OmegaToEuler(omegaEstRot(1:3:end), omegaEstRot(2:3:end), omegaEstRot(3:3:end));
[Model.lonEulerSig Model.latEulerSig Model.rateEulerSig ...
 Model.lonLatCorr Model.lonRateCorr Model.latRateCorr] = OmegaSigToEulerSig(omegaEstRot(1:3:end), omegaEstRot(2:3:end), omegaEstRot(3:3:end), Model.covarianceRot);

% Calculate rotation rates
Model.omegaX                                     = omegaEstRot(1:3:end);
Model.omegaY                                     = omegaEstRot(2:3:end);
Model.omegaZ                                     = omegaEstRot(3:3:end);

 
% Calculate rotation uncertainties
dEst                                             = sqrt(diag(Model.covarianceRot));
Model.omegaXSig                                  = dEst(1:3:end);
Model.omegaYSig                                  = dEst(2:3:end);
Model.omegaZSig                                  = dEst(3:3:end);

% Add zeros to the end of partials to account for the field rotation
% Calculate fault slip rates...
dEst                                             = Partials.slip * omegaEstRot(:);
Model.ssRate                                     = dEst(1:3:end);
Model.dsRate                                     = dEst(2:3:end);
Model.tsRate                                     = dEst(3:3:end);
% ...and uncertainties
dEst                                             = sqrt(diag(Partials.slip*Model.covarianceRot*Partials.slip'));
Model.ssRateSig                                  = dEst(1:3:end);
Model.dsRateSig                                  = dEst(2:3:end);
Model.tsRateSig                                  = dEst(3:3:end);

% Calculate forward components of velocity field
dEst                                             = (Partials.rotation-Partials.elastic * Partials.slip)*omegaEstRot;
Model.eastVel                                    = dEst(1:3:end);
Model.northVel                                   = dEst(2:3:end);
Model.upVel                                      = dEst(3:3:end);

% Rotational velocities
dEst                                             = Partials.rotation*omegaEstRot;
Model.eastRotVel                                 = dEst(1:3:end);
Model.northRotVel                                = dEst(2:3:end);
Model.upRotVel                                   = dEst(3:3:end);

% Velocities due to elastic strain accumulation on faults
dEst                                             = (Partials.elastic * Partials.slip)*omegaEstRot;
Model.eastDefVel                                 = dEst(1:3:end);
Model.northDefVel                                = dEst(2:3:end);
Model.upDefVel                                   = dEst(3:3:end);

% Calculate triangular slip rates...
[Model.trislipS,...
 Model.trislipD,...
 Model.trislipT,...
 Model.trislipSSig,...
 Model.trislipDSig,...
 Model.trislipTSig,...
 Model.trislipSBlock,...
 Model.trislipDBlock,...
 Model.trislipTBlock]                            = deal(zeros(size(Patches.v, 1), 1));
[Model.eastTriVel,...
 Model.northTriVel,... 
 Model.upTriVel]                                 = deal(zeros(size(Station.lon, 1), 1));     

if numel(omegaEstTriSlip) > 0
   if isempty(Command.triRake)
      Model.trislipS(triS)                       = omegaEstTriSlip(2*triS-1);
      Model.trislipD(triD)                       = omegaEstTriSlip(2*triD-0);
      Model.trislipT(triT)                       = omegaEstTriSlip(2*triT-0);
   else
      elrakes                                    = Command.triRake - rad2deg(Partials.tristrikes);
      Model.trislipS                             = omegaEstTriSlip.*cosd(elrakes);
      Model.trislipD                             = omegaEstTriSlip.*sind(elrakes);
   end
   % Calculate triangular slip uncertainties
   dEst                                          = sqrt(diag(Model.covarianceTriSlip));
   Model.trislipSSig(triS)                       = dEst(2*triS-1);
   Model.trislipDSig(triD)                       = dEst(2*triD-0);
   Model.trislipTSig(triT)                       = dEst(2*triT-0);
   
   % Calculate triangular slip rates as determined by block motion
%   dEst                                          = Partials.trislip * omegaEstRot(:);
%   Model.trislipSBlock                           = dEst(1:3:end);
%   Model.trislipDBlock                           = dEst(2:3:end);
%   Model.trislipTBlock                           = dEst(3:3:end);
   
   % Velocities due to elastic strain accumulation on triangular elements
   dEst                                          = Partials.tri*omegaEstTriSlip;
   Model.eastTriVel                              = dEst(1:3:end);
   Model.northTriVel                             = dEst(2:3:end);
   Model.upTriVel                                = dEst(3:3:end);
   
   % Modify model velocities
   Model.eastVel                                 = Model.eastVel - Model.eastTriVel;
   Model.northVel                                = Model.northVel - Model.northTriVel;
   Model.upVel                                   = Model.upVel - Model.upTriVel;
end

% Velocities due to internal block strains
if Command.strainMethod > 0
   dEst                                          = Partials.strain*omegaEstStrain;
   Model.eastStrainVel                           = dEst(1:2:end);
   Model.northStrainVel                          = dEst(2:2:end);
   Model.upStrainVel                             = zeros(size(Model.northStrainVel));
   
   % Modify model velocities
   Model.eastVel                                 = Model.eastVel + Model.eastStrainVel;
   Model.northVel                                = Model.northVel + Model.northStrainVel;
   Model.upVel                                   = Model.upVel + Model.upStrainVel;
else
   Model.eastStrainVel                           = zeros(size(Station.eastVel));
   Model.northStrainVel                          = zeros(size(Station.eastVel));
   Model.upStrainVel                             = zeros(size(Station.eastVel));
end   

% Calculate the residual velocites (done after consideration of strain components)
Model.eastResidVel                               = Station.eastVel - Model.eastVel;
Model.northResidVel                              = Station.northVel - Model.northVel;
Model.upResidVel                                 = Station.upVel - Model.upVel;

   
% Assign strain rates and uncertianties
Model.lonStrainSig                               = zeros(size(Model.omegaX));
Model.latStrainSig                               = zeros(size(Model.omegaX));
Model.eLonLon                                    = zeros(size(Model.omegaX));
Model.eLonLat                                    = zeros(size(Model.omegaX));
Model.eLatLat                                    = zeros(size(Model.omegaX));
Model.eLonLonSig                                 = zeros(size(Model.omegaX));
Model.eLonLatSig                                 = zeros(size(Model.omegaX));
Model.eLatLatSig                                 = zeros(size(Model.omegaX));
switch Command.strainMethod
   case 1 % 3, 3 parameter estimation 
      dEst                                       = omegaEstStrain;
      m1                                         = dEst(1:3:end);
      m2                                         = dEst(2:3:end);
      m3                                         = dEst(3:3:end);
      dEst                                       = sqrt(diag(Model.covarianceStrain));
      m1Sig                                      = dEst(1:3:end);
      m2Sig                                      = dEst(2:3:end);
      m3Sig                                      = dEst(3:3:end);
      strainBlockIdx                             = strainBlockIdx(3:3:end)/3;
      Model.eLonLon(strainBlockIdx)              = m1;
      Model.eLonLat(strainBlockIdx)              = m2;
      Model.eLatLat(strainBlockIdx)              = m3;
      Model.eLonLonSig(strainBlockIdx)           = m1Sig;
      Model.eLonLatSig(strainBlockIdx)           = m2Sig;
      Model.eLatLatSig(strainBlockIdx)           = m3Sig;
   case 2 % 5, 6 parameter estimation
      dEst                                       = omegaEstStrain;
      m1                                         = dEst(1:6:end);
      m2                                         = dEst(2:6:end);
      m3                                         = dEst(3:6:end);
      m4                                         = dEst(4:6:end);
      m5                                         = dEst(5:6:end);
      m6                                         = dEst(6:6:end);
      dEst                                       = sqrt(diag(Model.covarianceStrain));
      m1Sig                                      = dEst(1:6:end);
      m2Sig                                      = dEst(2:6:end);
      m3Sig                                      = dEst(3:6:end);
      m4Sig                                      = dEst(4:6:end);
      m5Sig                                      = dEst(5:6:end);
      m6Sig                                      = dEst(6:6:end);
      strainBlockIdx                             = strainBlockIdx(6:6:end)/6;
%     Model.lonStrain(strainBlockIdx)            = 
%     Model.latStrain(strainBlockIdx)            = 
      Model.eLonLon(strainBlockIdx)              = m1.*m2./m4;
      Model.eLonLat(strainBlockIdx)              = m2;
      Model.eLatLat(strainBlockIdx)              = m5;
      Model.eLonLonSig(strainBlockIdx)           = m1Sig.*m2Sig./m4Sig;
      Model.eLonLatSig(strainBlockIdx)           = m2Sig;
      Model.eLatLatSig(strainBlockIdx)           = m5Sig;
   case 3 % 3, 4 parameter estimation
      dEst                                       = omegaEstStrain;
      m1                                         = dEst(1:4:end);
      m2                                         = dEst(2:4:end);
      m3                                         = dEst(3:4:end);
      m4                                         = dEst(4:4:end);
      dEst                                       = sqrt(diag(Model.covarianceStrain));
      m1Sig                                      = dEst(1:4:end);
      m2Sig                                      = dEst(2:4:end);
      m3Sig                                      = dEst(3:4:end);
      m4Sig                                      = dEst(4:4:end);
      strainBlockIdx                             = strainBlockIdx(4:4:end)/4;
      Model.eLonLon(strainBlockIdx)              = m1.*m2./m3;
      Model.eLonLat(strainBlockIdx)              = m2;
      Model.eLatLat(strainBlockIdx)              = m4;
      Model.eLonLonSig(strainBlockIdx)           = m1Sig.*m2Sig./m3Sig;
      Model.eLonLatSig(strainBlockIdx)           = m2Sig;
      Model.eLatLatSig(strainBlockIdx)           = m4Sig;
   case 4 % directed forward search method
   
end

fprintf('done.\n')

if Command.nIter > 1
   fprintf('Carrying out Monte Carlo simulation of slip rates and uncertainties...')
   warning off all
   nIter = Command.nIter;
   
   [ssrate, dsrate, tsrate] = deal(zeros(size(Partials.slip, 1)/3, nIter));
   
   for i = 1:nIter
      fprintf('\n  Working on iteration %d of %d...', i, nIter)
      sigs = reshape([Station.eastSig Station.northSig]', 2*numel(Station.lon), 1);
      newSigs = sigs.*(2*rand(2*numel(Station.lon), 1) - 1);
      d(1:2*numel(Station.lon)) = d(1:2*numel(Station.lon)) + newSigs;
   
      %Command.triCons = [0 150];
      Command.triCons = [];
      [R, W, d]                                        = deal(full(R), full(W), full(d));
      Model.covariance                                 = inv(R'*W*R);
      if isempty(Command.triCons)
         omegaEst                                      = Model.covariance*R'*W*d;
      else
         [minconstr, maxconstr]                        = deal(-inf(size(R, 2), 1), inf(size(R, 2), 1));
         minconstr(szslip(2)+1:szslip(2)+sztri(2))     = Command.triCons(1);
         maxconstr(szslip(2)+1:szslip(2)+sztri(2))     = Command.triCons(2);
         [omegaEst, resnorm, res, exf, output, lambda] = lsqlin(R, d, [], [], [], [], minconstr, maxconstr);
      end
      
      % extract different parts of the state vector
      omegaEstRot                                      = omegaEst(1:size(Partials.slip, 2)); % rotation parameters
      
      
      % Add zeros to the end of partials to account for the field rotation
      % Calculate fault slip rates...
      dEst                                             = Partials.slip * omegaEstRot(:);
      ssrate(:, i)                                  = dEst(1:3:end); 
      dsrate(:, i)                                  = dEst(2:3:end); 
      tsrate(:, i)                                  = dEst(3:3:end);
   end
   
   Model.ssRate                                     = mean(ssrate, 2);
   Model.dsRate                                     = mean(dsrate, 2);
   Model.tsRate                                     = mean(tsrate, 2);
   
   Model.ssRateSig                                  = std(ssrate, 0, 2);
   Model.dsRateSig                                  = std(dsrate, 0, 2);
   Model.tsRateSig                                  = std(tsrate, 0, 2);

   fprintf('done.\n')
end

fprintf('Writing output...')
% Write output
WriteOutput(Segment, Patches, Station, Block, Command, Model, Partials.tristrikes);
%save BlocksStore.mat
system(sprintf('cp %s .%s%s%s.', commandFile, filesep, runName, filesep));
if exist('tempkernels.mat') == 2;
   system(sprintf('mv tempkernels.mat .%s%s%skernels.mat', filesep, runName, filesep));
end
save([runName filesep 'ModelStruct.mat'], 'Model')
fprintf('done.  All files saved to .%s%s.\n', filesep, runName)
