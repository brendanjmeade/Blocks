function Model = MogiResults(Model);

Model.mogiDeltaV = 1e-3*Model.omegaEstMogi; % Converting to cubic meters/yr
Model.mogiDeltaVSig = 1e-3*sqrt(diag(Model.covarianceMogi)); % Converting to cubic meters/yr