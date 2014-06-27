function G = GetElasticPartials(F, S)
% GetElasticPartials   Calculates elastic displacement partial derivatives.
%    G = GetElasticPartials(F, S) calculates the elastic displacement partial
%    derivatives based on the Okada formulation, using the source and receiver
%    geometries defined in structures F and S.  F gives the fault geometry, in
%    Blocks' .segment format, and S gives the station geometry, in Blocks' 
%    .sta.data format.  Before calculating the partials for each segment, a local
%    oblique Mercator project is done.
%

% Hard-coded Poisson ratio
nu                                          = 0.25;

% Count input features
nSta                                        = numel(S.lon);
nSeg                                        = numel(F.lon1);

% Allocate space for partials
G                                           = zeros(3*nSta, 3*nSeg);
[v1 v2 v3]                                  = deal(cell(1, nSeg));

% Calculate projected Cartesian coordinates
[F, S]                                      = ProjectSegCoords(F, S);

% Calculate fault parameters in Okada form
[strike, L, W, ofx, ofy, ofxe, ofye, ...
               tfx, tfy, tfxe, tfye]        = fault_params_to_okada_form(F.px1, F.py1, F.px2, F.py2, deg_to_rad(F.dip), F.lDep, F.bDep);


% Loop through each segment and calculate displacement
parfor (i = 1:nSeg)

    % Displacements due to unit slip components
    [ves vns vus...
     ved vnd vud...
     vet vnt vut]                           = okada_partials(ofx(i), ofy(i), strike(i), F.lDep(i), deg_to_rad(F.dip(i)), L(i), W(i), 1, 1, 1, S.fpx(:, i), S.fpy(:, i), nu);
     
    v1{i}                                   = reshape([ves vns vus]', 3*nSta, 1);
    v2{i}                                   = reshape(sign(90 - F.dip(i)).*[ved vnd vud]', 3*nSta, 1);
    v3{i}                                   = reshape((F.dip(i) - 90 == 0).*[vet vnt vut]', 3*nSta, 1); 

end

% Place cell arrays into the partials matrix
G(:, 1:3:end)                               = cell2mat(v1);
G(:, 2:3:end)                               = cell2mat(v2);
G(:, 3:3:end)                               = cell2mat(v3);

% Rotate displacement vectors to project back into global coordinate system
% (removing effects of local oblique Mercator projection)
G                                           = xyz2enumat(G, -F.strike + 90);
