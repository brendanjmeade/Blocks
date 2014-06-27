function Blocks_noise(commandFile)

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
   [apLons, apLats, apRates]                     = deal(Block.eulerLon(blockConIdx), Block.eulerLat(blockConIdx), Block.rotationRate(blockConIdx).*1e6);
   [apbx, apby, apbz]                            = EulerToOmega(apLons, apLats, apRates);
   blockCons(1:3:end)                            = apbx;
   blockCons(2:3:end)                            = apby;
   blockCons(3:3:end)                            = apbz;
   [apbsx, apbsy, apbsz]                         = epoles_cov_to_omega_cov(apbx, apby, apbz, diag([deg_to_rad(Block.eulerLatSig(blockConIdx)), deg_to_rad(Block.eulerLonSig(blockConIdx)), deg_to_rad(Block.rotationRateSig(blockConIdx))]));
   blockConSigs(1:3:end)                         = apbsx;
   blockConSigs(2:3:end)                         = apbsy;
   blockConSigs(3:3:end)                         = apbsz;
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
   [Wtri, Wseg, beta]                            = SmoothEdges(share, Patches, Segment, Partials.slip, Command);                 
   Wtri                                          = Wtri(:, colkeep); % eliminate the entries corresponding to empty partials
   Wtri                                          = Wtri(colkeep, :); % eliminate the entries corresponding to empty partials
   Wseg(3:3:end, :)                              = [];
   beta(3:3:end)                                 = [];
   Wconst                                        = Wconst(:, colkeep);
   Wcons                                         = [Wconss Wconst];
   Wcons(3:3:end, :)                             = [];
else
   [Wseg, Wtri, Wcons, beta, Ztri]               = deal([]);
end
sztri                                            = size(Partials.tri);

fprintf('\n  Assembling the Jacobian...')
% Assemble Jacobian

% First determine submatrix sizes and indices

% Submatrix row counts
rsta                                             = 2*numel(Station.lon); % station rows
rbcons                                           = size(Partials.blockCon, 1); % Block constraint rows
rscons                                           = size(Partials.slipCon, 1); % Slip constraint rows
rtriw                                            = sztri(2); % Triangular smoothing rows
rtriz                                            = size(Ztri, 1); % Triangular edge constraint rows

% Submatrix row indices
rnum                                             = [rsta rbcons rscons rtriw rtriz];
ridx                                             = cumsum([0 rnum]);

% Calculate noise partials as a big identity matrix
Partials.noise                                   = eye(rsta);
%Partials.noise                                   = eye(sum(rnum));

% Submatrix column counts
cblock                                           = szslip(2); % Block/slip columns
ctri                                             = rtriw; % Triangle columns
cnoise                                           = size(Partials.noise, 2); % Noise columns

% Submatrix column indices
cnum                                             = [cblock ctri cnoise];
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
R(rows{1, 3}, cols{1, 3})                        = Partials.noise;
%R(:, cols{1, 3})                                 = Partials.noise;
R(rows{2, 1}, cols{2, 1})                        = Partials.blockCon;
R(rows{3, 1}, cols{3, 1})                        = Partials.slipCon;
R(rows{4, 2}, cols{4, 2})                        = Wtri;
R(rows{5, 2}, cols{5, 2})                        = Ztri;

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

wd                                               = zeros(3*numel(Station.lon), 1);
wd(1:3:end)                                      = 1./Station.eastSig.^2;
wd(2:3:end)                                      = 1./Station.northSig.^2;
wd(3:3:end)                                      = 1./Station.upSig.^2;
wd                                               = wd(:);
wd                                               = 1*wd(rowkeep);
wd                                               = [wd ; Command.blockConWgt.*1./blockConSigs.^2]; % Add block motion constraints
wd                                               = [wd ; Command.slipConWgt.*1./(slipSigs(slipConIdx)).^2]; % Add slip rate constraints
wd                                               = [wd ; beta]; % add the triangular smoothing vector
wd                                               = [wd ; ones(size(Wcons, 1), 1)]; % add the triangular kinematic consistency vector
wd                                               = [wd ; 1e5*ones(size(Ztri, 1), 1)]; % add the zero edge vector
Wd                                               = spdiags(1./wd, 0, numel(wd), numel(wd)); % assemble into a matrix
% Make a copy for iteration
wdi                                              = wd;

slipSig                                          = 1;
noiseSig                                         = 1;
sigRatio                                         = slipSig./noiseSig;
wm                                               = ones(cidx(end), 1);
wm(cols{1, 3}(1:2:end))                          = noiseSig*Station.eastSig.^2;
wm(cols{1, 3}(2:2:end))                          = noiseSig*Station.northSig.^2;
%w2(rsta+(1:sum(slipToggle)))                     = noiseSig*Command.slipConWgt*slipSigs(slipConIdx);
Wm                                               = diag(wm);

slipRange = 10.^[-3:0.1:4];
slipRange = 10.^[-9:1:2];
slipRange = 10.^-6;
%noiseRange = 10.^[-6:.25:2];
%[L, mlik] = maxlike2d(slipRange, noiseRange, cblock + ctri, R, Wd, Wm, d);
%[~, midx] = min(-L);
%slipSig = slipRange(midx);
%keyboard
noiseRange = 10.^[-3:1:5];
dvec = zeros(2*numel(Station.lon), 1);
dvec(1:2:end) = Station.eastSig;
dvec(2:2:end) = Station.northSig; 
dvec = [dvec];% blockConSigs(:); slipSigs(slipConIdx); beta; ones(size(Wcons, 1), 1); ones(size(Ztri, 1), 1)];
dvec = dvec + 0.01*randn(size(dvec));
[slipSig, noiseSig] = optparams(slipRange, noiseRange, cblock + ctri, R, Wd, Wm, d, dvec);
%slipSig = 1e-3;
%noiseSig = 1;
%wm                                               = slipSig*ones(cidx(end), 1);
%wm(cols{1, 3})                                   = noiseSig*ones(length(cols{1, 3}), 1);
wm(cols{1, 1})                                   = slipSig*wm(cols{1, 1});
wm(cols{1, 3})                                   = noiseSig*wm(cols{1, 3});
Wm                                               = diag(wm);

% Make a copy for iteration
wmi                                              = wm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try distance weighted correlation-based covariance matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%w2                                               = eye(cidx(end));
%disweight = 5;
%for i = 1:(rsta/2)
%   di = distance(Station.lat(i), Station.lon(i), Station.lat, Station.lon);
%   dw = 1./(1+di).^disweight;
%   w2(2*i-1, 1:2:rsta) = dw;
%   w2(2*i-0, 2:2:rsta) = dw;
%end
%W2 = noiseSig*w2;
%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('done.\n')

% Estimate rotation vectors, Euler poles and slip rates
fprintf('Doing the inversion...\n')

% If you receive out of memory errors during the inversion, comment out the next line.
% This will likely allow the inversion to be carried out, but it will be single-threaded
[R, Wd, d]                                        = deal(full(R), full(Wd), full(d));

%Model.covariance                                 = (R*W*R')\eye(size(R, 1));
%omegaEst                                         = W*R'*Model.covariance*d;
%Model.covariance                                 = (W + R*R')\eye(size(R, 1));
%omegaEst                                         = R'*Model.covariance*d;

inv1                                             = (Wd + R*Wm*R')\eye(size(R, 1));
omegaEst                                         = Wm*R'*inv1*d;
res                                              = d - R*omegaEst;
%mlik                                             = 1./length(res)*res'*inv1*res;
%logdet                                           = sum(log(diag(chol(Wd + R*Wm*R'))));
%loglik                                           = -0.5*(length(res) - length(res)*log(length(res))) - 0.5*logdet - 0.5*length(res)*log(res'*inv1*res)
%toler                                            = abs(mlik - 1);
%
%% Set up recursion for maximum likelihood scaling of model covariance
%niter = 1;
%maxiter = 100;
%while toler > 1e-10 & niter <= maxiter
%   noiseSig                                      = noiseSig*mlik;
%   slipSig                                       = slipSig/mlik;
%   wm(cols{1, 1})                                = sigRatio*slipSig;
%   wm(cols{1, 3})                                = noiseSig*wmi(cols{1, 3});
%   Wm                                            = diag(wm);
%   wd(rows{1, 1})                                = noiseSig*wdi(rows{1, 1});
%   Wd                                            = diag(wd);
%   inv1                                          = (Wd + R*Wm*R')\eye(size(R, 1));
%   omegaEst                                      = Wm*R'*inv1*d;
%   res                                           = d - R*omegaEst;
%   mlik                                          = 1./length(res)*res'*inv1*res
%   logdet                                        = sum(log(diag(chol(Wd + R*Wm*R'))));
%   loglik                                        = -0.5*(length(res) - length(res)*log(length(res))) - 0.5*logdet - 0.5*length(res)*log(res'*inv1*res)
%   toler                                         = abs(mlik - 1);
%   niter                                         = niter + 1;
%end
%
%keyboard
Model.covariance                                 = Wm*R'*inv1*dot(res, res)*inv1*R*Wm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try different methods of inversion for underdetermined system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pinv1                                            = pinv(W + R*W2*R');
%omegaEst1                                         = W2*R'*pinv1*d;
%res1                                              = d - R*omegaEst1;
%norm(omegaEst1)
%
%[RQ, RR]                                          = qr((W + R*W2*R')');
%omegaEst2                                         = W2*R'*RQ*(RR'\d);
%res2                                              = d - R*omegaEst2;
%norm(omegaEst2)
%
%Model.covariance                                 = W2*R'*inv1*dot(res, res)*inv1*R*W2;
%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nInversion is done.\n')

fprintf('Calculating model results...')
% extract different parts of the state vector
omegaEstRot                                      = omegaEst(1:szslip(2)); % rotation parameters
omegaEstTriSlip                                  = omegaEst(szslip(2) + (1:sztri(2))); % triangular slips
omegaEstNoise                                    = omegaEst(szslip(2) + sztri(2) + 1:end); % noise parameters

% extract different parts of the model covariance matrix
Model.covarianceRot                              = Model.covariance(1:szslip(2), 1:szslip(2));
Model.covarianceTriSlip                          = Model.covariance(szslip(2) + (1:sztri(2)), szslip(2) + (1:sztri(2)));
Model.covarianceNoise                            = Model.covariance(szslip(2) + sztri(2) + 1:end, szslip(2) + sztri(2) + 1:end);

% Calculate Euler poles and uncertainties
[Model.rateEuler,...
 Model.lonEuler,...
 Model.latEuler]                                 = OmegaToEuler(omegaEstRot(1:3:end), omegaEstRot(2:3:end), omegaEstRot(3:3:end));
[Model.lonEulerSig, Model.latEulerSig,...
 Model.rateEulerSig, Model.lonLatCorr,...
 Model.lonRateCorr, Model.latRateCorr]           = OmegaSigToEulerSig(omegaEstRot(1:3:end), omegaEstRot(2:3:end), omegaEstRot(3:3:end), Model.covarianceRot);

% Calculate rotation rates
Model.omegaX                                     = omegaEstRot(1:3:end);
Model.omegaY                                     = omegaEstRot(2:3:end);
Model.omegaZ                                     = omegaEstRot(3:3:end);
 
% Calculate rotation uncertainties
dEst                                             = sqrt(diag(Model.covarianceRot));
Model.omegaXSig                                  = dEst(1:3:end);
Model.omegaYSig                                  = dEst(2:3:end);
Model.omegaZSig                                  = dEst(3:3:end);
%
%% ****** FAULT SLIP RATES USING MONTE CARLO ******
%niter = 100;
%omegaEstMC = zeros(size(omegaEst, 1), niter);
%slipsMC = zeros(3*numel(Segment.lon1), niter);
%sigs = reshape([Station.eastSig Station.northSig]', 2*numel(Station.lon), 1);
%dnew = d;
%noiseFac = 2;
%for i = 1:niter
%   newSigs = sigs.*(2*noiseFac*rand(2*numel(Station.lon), 1) - noiseFac);
%   dnew(1:2*numel(Station.lon)) = d(1:2*numel(Station.lon)) + newSigs;
%   omegaEstMC(:, i) = W2*R'*inv1*dnew;
%   omegaEstRotMC = omegaEstMC(1:szslip(2), i);
%   slipsMC(:, i) = Partials.slip * omegaEstRotMC;
%end
%slipsmean = mean(slipsMC, 2);
%dEst = slipsmean;
%Model.ssRate                                     = dEst(1:3:end);
%Model.dsRate                                     = dEst(2:3:end);
%Model.tsRate                                     = dEst(3:3:end);
%slipssig  = std(slipsMC, 0, 2);
%dEst = slipssig;
%Model.ssRateSig                                  = dEst(1:3:end);
%Model.dsRateSig                                  = dEst(2:3:end);
%Model.tsRateSig                                  = dEst(3:3:end);


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

% Extract the estimated noise as velocity components
Model.eastStrainVel                              = omegaEstNoise(1:2:rsta);
Model.northStrainVel                             = omegaEstNoise(2:2:rsta);

% And correct modeled velocities
Model.eastVel                                    = Model.eastVel + Model.eastStrainVel;
Model.northVel                                   = Model.northVel + Model.northStrainVel;

% Calculate the residual velocites (done after consideration of strain components)
Model.eastResidVel                               = Station.eastVel - Model.eastVel;
Model.northResidVel                              = Station.northVel - Model.northVel;
Model.upResidVel                                 = Station.upVel - Model.upVel;
   
fprintf('done.\n')

fprintf('Writing output...')
% Write output
WriteOutputNoise(Segment, Patches, Station, Block, Command, Model, Partials.tristrikes);
% Uncomment the next line to save a giant file with all variables
%save BlocksStore.mat
system(sprintf('cp %s .%s%s%s.', commandFile, filesep, runName, filesep));
% Move the elastic kernels to the output directory
if exist('tempkernels.mat') == 2;
   system(sprintf('mv tempkernels.mat .%s%s%skernels.mat', filesep, runName, filesep));
end
fprintf('done.  All files saved to .%s%s.\n', filesep, runName)
