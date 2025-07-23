function [node, elem] = loadMesh
%% GMSH Generate mesh for surface
%   [NODE, ELEM] = GMSH returns mesh information in the data
%   structure in NODE and ELEM
%
%   Copyright (C) Hailong Guo
%   11/17/2016

MeshDirectory = 'mesh/' ;

%% Load mesh
NodeName = [MeshDirectory 'node'];
mesh = load(NodeName);
node = mesh.node;

ElemName = [MeshDirectory 'elem'];
mesh = load(ElemName);
elem = mesh.elem;

