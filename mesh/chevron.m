function mesh = chevron(nx, ny)
%CHEVRON Generate chevron pattern uniform mesh
%   [NODE, ELEM] = CHEVRON(NX, UY) generate chevron pattern uniform mesh
%   with mesh size in x-direction 1/nx and mesh size in y-direction 1/ny
%
%   Copyright (C) Hailong Guo
%   03/13/2014

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

%% Generate element information
elem = zeros(2*nx*ny, 3);
for i = 1:ny
    for j = 1:4:2*nx
        elem((i-1)*2*nx+j, :) = [1, nx+3, nx+2,]+(i-1)*(nx+1)+floor((j-1)/2);
    end
    for j = 2:4:2*nx
        elem((i-1)*2*nx+j, :) = [1, 2, nx+3,]+(i-1)*(nx+1)+floor((j-1)/2);
    end
    for j = 3:4:2*nx
        elem((i-1)*2*nx+j, :) = [2, 3, nx+3,]+(i-1)*(nx+1)+floor((j-3)/2);
    end
    for j = 4:4:2*nx
        elem((i-1)*2*nx+j, :) = [3, nx+4, nx+3,]+(i-1)*(nx+1)+floor((j-3)/2);
    end
end

%% contruct the data structure
mesh = struct('node', node, 'elem', elem);
