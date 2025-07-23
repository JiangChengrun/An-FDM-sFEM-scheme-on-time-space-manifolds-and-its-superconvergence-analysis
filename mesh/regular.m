function mesh = regular(nx, ny, rect)
%% REGULAR Generate uniform mesh of regular pattern
%   mesh = REGULAR(NX, NY) returns mesh
%
%   Copyright (C) Hailong Guo
%   Date: 06/01/2017

%% Initialization
if nargin == 1
    ny = nx;
    x0 = 0; x1 = 1;
    y0 = 0; y1 = 1;
elseif nargin == 2
    x0 = 0; x1 = 1;
    y0 = 0; y1 = 1;
elseif nargin == 3
    x0 = rect(1); x1 = rect(2);
    y0 = rect(3); y1 = rect(4);
end

%% Generate vertice coordinate
hx = (x1-x0)/nx;
hy = (y1-y0)/ny;
[x,y] = meshgrid(x0:hx:x1,y0:hy:y1);
node(:,1) = x(:);
node(:,2) = y(:);
node(:,3) = 0;

%% Generate element numbering
elem = zeros(2*nx*ny, 3);
for i = 1:nx
    for j = 1:ny
        elem((i-1)*2*ny+2*j-1,:) = [1 ny+2 ny+3]+(i-1)*(ny+1)+j-1;
        elem((i-1)*2*ny+2*j,:) = [1 ny+3 2]+(i-1)*(ny+1)+j-1;
    end
end

%% Main data structure
mesh = struct('node',node,'elem',elem, 'nx', nx, 'ny', ny);
