function WriteMogi(mogifile, Mogi, Model)
% WRITEMOGI  Writes Mogi sources
%    WRITEMOGI(file) writes a 5-column array to the 
%    specified file, of the format:
%
%    longitude | latitude | depth (km) | \DeltaV (km^3) | \DeltaVSigma (km^3)
%
%    The file is given a one-line header, and 
%    the columns are space-separated.
%


fid        = fopen(mogifile, 'w');
fprintf(fid, 'Longitude Latitude Depth DV_Flag DV DVS\r');
fprintf(fid, '%g %g %g %g %g %g\r', [Mogi.lon, Mogi.lat, Mogi.dep, ones(size(Mogi.lon)), Model.mogiDeltaV, Model.mogiDeltaVSig]');
fclose(fid);