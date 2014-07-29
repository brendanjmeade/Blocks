function [p, pp] = readmshp(file, pp)
% readmshp  Reads a mesh property file. 
%   [P, PP] = readmshp(FILE) reads mesh names and properties from the 
%   specified FILE. The .mshp file contains groups of 3 lines and
%   specifies the file name for a mesh along with its properties:
%
%   1. Full path to mesh file (.msh or .mat), readable using ReadPatches
%   2. Smoothing coefficient (double)
%   3. A priori edge constraints as a triplet, with values corresponding
%      to conditions on the [updip downdip lateral] edges of the mesh. 
%      Numerical values can be 0 (no constraint), 1 (creeping), or 2 
%      (fully coupled).
%   4. Full path to associated a priori slip file (.mat) (blank if none)
%   
%   File contents are returned to structure P, containing the actual 
%   triangular meshes, and PP, containing properties of the meshes. 
%
%   [P, COMMAND] = readmshp(FILE, COMMAND) updates an input COMMAND structure
%   with mesh properties. 
%

% Read file contents
fid = fopen(file, 'r');
c = textscan(fid, '%s\n%f\n%f%f%f\n%s\n');
pp.patchFileNames = char(c{1});

% Read patches
p = ReadPatches(pp.patchFileNames);

% Define mesh properties
if exist('command', 'var')
   pp = command;
end
pp.triSmooth = c{2}(:)';
pp.triEdge = reshape([c{3} c{4} c{5}]', 1, 3*size(c{3}, 1));
pp.slipFileNames = char(c{6});