function Patches                 = ReadPatches(filenames);
% 
% ReadPatches.m
%
% This function reads in any triangulated patch files specified by
% Segment.other1 and returns their coordinates to Patches.
%
% Arguments:
%   filenames    :  string containing all filenames of patch files
%
% Returned variables:
%   Patches      :  a structure containing:
%          .c    : all vertex coordinates (n x 3)
%          .v    : element vertex indices (m x 3)
%          .nEl  : number of elements in each patch file
%          .nc   : number of coordinates in each patch file
%
Patches.c                        = [];
Patches.v                        = [];
Patches.nEl                      = [];
Patches.nc                       = [];

if numel(filenames) > 0
   if size(filenames, 1) == 1
      spaces                     = [0 findstr(filenames, ' ') length(filenames)+1];
      nfiles                     = length(spaces) - 1;
   else
      nfiles                     = size(filenames, 1);
   end
   
   for i = 1:nfiles
      if exist('spaces', 'var')
         filename                = filenames(spaces(i)+1:spaces(i+1)-1);
      else
         filename                = strtrim(filenames(i, :));
      end
      
      if filename(end-3:end) == '.msh'
         [c, v]                  = msh2coords(filename);
      elseif filename(end-3:end) == '.mat'
         load(filename, 'c', 'v')
      end
      crossd                     = cross([c(v(:, 2), :) - c(v(:, 1), :)], [c(v(:, 3), :) - c(v(:, 1), :)]);
      negp                       = find(crossd(:, 3) < 0);
      [v(negp, 2), v(negp, 3)]   = swap(v(negp, 2), v(negp, 3));
      Patches.v                  = [Patches.v; v + sum(size(Patches.c, 1))];
      Patches.c                  = [Patches.c; c];
      Patches.nEl                = [Patches.nEl; size(v, 1)];
      Patches.nc                 = [Patches.nc; size(c, 1)];
   end
end