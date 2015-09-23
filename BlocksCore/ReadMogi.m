function Mogi = ReadMogi(mogifile)
% READMOGI  Reads in coordinates of Mogi sources
%    M = READMOGI(file) reads a 6-column array from the 
%    specified file, of the format:
%
%    longitude, latitude, depth (km), toggle, \DeltaV (m^3), \DeltaVSigma (m^3)
%
%    The file is assumed to have a one-line header, and 
%    the columns are assumed to be comma-separated.
%

if ~isempty(mogifile)
   fid        = fopen(mogifile, 'r');
   m          = textscan(fid, '%f, %f, %f, %f, %f, %f\n', 'headerlines', 1);
   Mogi.lon   = m{1};
   Mogi.lat   = m{2};
   Mogi.dep   = m{3};
   Mogi.dvtog = m{4};
   Mogi.dv    = m{5};
   Mogi.dvSig = m{6};
else
   Mogi.lon   = [];
   Mogi.lat   = [];
   Mogi.dep   = [];
   Mogi.dvtog = [];
   Mogi.dv    = [];
   Mogi.dvSig = [];
end   