function [node, elem] = EnzensbergerSternMesh
%% EnzensbergerSternMesh Generate mesh using EnzensbergersternMesh 
%
%   Copyright (C) Hailong Guo
%   11/27/2016

%% Addpath 
addpath('~/Src/matlab/distmesh');

%% Preparation
surfacedata = EnzensbergerSternSurface();
fd = @(p) 400*(p(:,1).^2.*p(:,2).^2 + p(:,2).^2.*p(:,3).^2 + p(:,1).^2.*p(:,3).^2) ...
- (1-p(:,1).^2-p(:,2).^2-p(:,3).^2).^3 - 40;

%% Generate data
figure(1)
[node,elem]=distmeshsurface(fd,@huniform,0.2,[-2.5,-2.5,-2.5; 2.5,2.5,2.5]);
node  = surfacedata.project(node);

%% Save data for reusing
save('mesh/node', 'node');
save('mesh/elem', 'elem');