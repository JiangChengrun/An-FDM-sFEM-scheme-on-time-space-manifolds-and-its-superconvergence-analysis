function mesh = readMesh(filename)
%% READMESH Read mesh generate by CGAL 
%   [NODE, ELEM] = READMESH(FILENMAE) read mesh genrated by CGAL and save the 
%   corrdinates of verticesin NODE and vertice connection information in ELEM
%
%   Copyright (C) Hailong Guo
%   $12/06/2016

%% Open the file
fid = fopen(filename, 'r');

%% Parse the input file
fgetl(fid);
tline = fscanf(fid, '%d %d %d', [1, 3]);
nbv = tline(1);
nbt = tline(2);
fgetl(fid);
node = fscanf(fid, '%f %f \n', [ 3, nbv])';
elem = fscanf(fid, '%d %d %d \n', [4, nbt])';
elem = elem(:, 2:4) + 1;

%% Construct the mesh structure
mesh.node = node;
mesh.elem = elem;