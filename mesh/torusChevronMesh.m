function mesh = torusChevronMesh(Nu, Nv, R, r)
%% TORUSUNIFORMMESH Generate chevron pattern mesh on torus
%   [NODE,ELEM] = TORUSUNIFORMMESH(NU, NV, R, R) return the mesh data of
%   chevron pattern mesh on torus
%
%   Copyright (C) Hailong Guo
%   11/24/2016

%% Initial Data

% Nu = 4;
% Nv = 4;
% R = 4;
% r = 1;

%% Initialize
u=(0:2*pi/Nu:2*pi)';
v=(0:2*pi/Nv:2*pi)';

%% Generate element connect information
elem=zeros(2*Nu*Nv,3);
count =0;
for i=1:Nu
    for j=1:Nv
        if (i<Nu && j<Nv)
            if mod(i,2) == 1
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j,i*Nv+j+1, (i-1)*Nv+j+1];
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j, i*Nv+j, i*Nv+j+1];
            else
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j,i*Nv+j, (i-1)*Nv+j+1];
                count = count + 1;
                elem(count,:)=[i*Nv+j, i*Nv+j+1, (i-1)*Nv+j+1];
            end
        elseif (i==Nu &&j<Nv)
            if mod(i,2) == 1
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j, j+1, (i-1)*Nv+j+1];
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j, j, j+1];
            else
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j, j, (i-1)*Nv+j+1];
                count = count + 1;
                elem(count,:)=[j, j+1, (i-1)*Nv+j+1];
            end
        elseif (i<Nu && j==Nv)
            if mod(i,2) == 1
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j,i*Nv+1, (i-1)*Nv+1];
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j, i*Nv+j, i*Nv+1];
            else
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j,i*Nv+j, (i-1)*Nv+1];
                count = count + 1;
                elem(count,:)=[i*Nv+j, i*Nv+1, (i-1)*Nv+1];
            end
        else
            if mod(i,2) == 1
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j, 1, (i-1)*Nv+1];
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j, j, 1];
            else
                count = count + 1;
                elem(count,:)=[(i-1)*Nv+j, j, (i-1)*Nv+1];
                count = count + 1;
                elem(count,:)=[j, 1, (i-1)*Nv+1];
            end
        end
    end
end

%% Generate node coordinate
node=zeros(Nu*Nv,3);
count = 0;
for i=1:Nu
    for j=1:Nv
        count = count + 1;
        node(count,:)=[(R+r*cos(v(j)))*cos(u(i)),(R+r*cos(v(j)))*sin(u(i)),r*sin(v(j))];
        %node(count,:) = [u(i),v(j),0];
    end
end

%% Construct the data structure
mesh = struct('node', node, 'elem', elem);