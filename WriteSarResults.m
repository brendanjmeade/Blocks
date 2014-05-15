function out = WriteSarResults(Sar, Model, runName)

out = [Sar.lon, Sar.lat, Sar.data, Model.Sar, Model.SarResid, Model.SarDef, Model.SarRot, Model.SarRamp, Model.SarShift, Model.SarTri, Model.SarStrain];
save([runName filesep 'Sar.pred'], 'out', '-ascii');
