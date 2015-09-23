function WriteMogi(mogifile, Mogi, Model)
% WRITEMOGI  Writes Mogi sources
%    WRITEMOGI(file) writes a 6-column array to the 
%    specified file, of the format:
%
%    longitude, latitude, depth (km), toggle, \DeltaV (m^3), \DeltaVSigma (m^3)
%
%    The file is given a one-line header, and 
%    the columns are comma-separated.
%


fid        = fopen(mogifile, 'w');
fprintf(fid, 'Longitude, Latitude, Depth, DV_Flag, DV, DVS\r');
fprintf(fid, '%g, %g, %g, %g, %g, %g\r', [Mogi.lon, Mogi.lat, Mogi.dep, ones(size(Mogi.lon)), Model.mogiDeltaV, Model.mogiDeltaVSig]');
fclose(fid);