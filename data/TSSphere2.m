function surfacedata = TSSphere2
%% ELLIOTTSURFACE  Data for Laplace-beltrami equation
%   SURFACEDATA = ELLIOTTSURFACE returns the exact solution and right hand
%   side function of the following equation
%         - \Delta_{\Gamma} u + u = f
%   on the surface given by the level set function phi. It is the exmaple
%   3 in the following refernece:
%
%   [1]. Chernyshenko, Alexey Y. ;  Olshanskii, Maxim A.  An adaptive
%        octree finite element method for PDEs posed on surfaces.
%         Comput. Methods Appl. Mech. Engrg.  291  (2015), 146--172.
%
%   Copyright (C) Guozhi Dong & Hailong Guo
%   02/21/2017

%% Main data structure
surfacedata = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'initmesh',@initmesh,...
    'meancurvature', @meancurvature, 'exactu', @exactu, 'f', @f, 'Du',@Du, ...
    'g_N0',@g_N0, 'g_N1', @g_N1,'projectNT', @projectNT);

%% Subfaction
% Inital mesh
    function [node, elem] = initmesh
        node = [];
        elem = [];
        warning('That function is not implemented');
    end

% level set function
    function s = phi(p)
        % level set function
        x = p(:,1); y = p(:,2); z = p(:,3);
        s = sqrt(x.^2+y.^2+z.^2)-1;
    end

% unit out normal vector
    function un = unitoutnormal(p)
        n = gradient(p);
        nn = sqrt(sum(n.^2,2));
        un = n./repmat(nn,1,3);
    end

% gradient of level set function
    function n = gradient(p)
        L = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
        n = p ./(L*ones(1,3));
    end

% project point on to the surface, this from is part of iFEM by Long Chen
    function node = project(p)
        % projection function
        d = phi(p);
        node = p - d*ones(1,3).*gradient(p);
    end


% right hand side function
    function s = f(p, ~)
        p = project(p);
        x = p(:,1); y = p(:,2); z = p(:,3);
        r = sqrt(x.^2+y.^2+z.^2);
        phix = x./r;
        phiy = y./r;
        %phiz = z./r;
        phixx = 1./r - x.^2./r.^3;
        phiyy = 1./r - y.^2./r.^3;
        phizz = 1./r - z.^2./r.^3;
        dn = phixx + phiyy + phizz;
        s = (2*phix.*phiy + dn.*(y.*phix + x.*phiy)) + x.*y;
    end

% exact solution
    function rhs = exactu(p, ~)
        p = project(p);
        x = p(:,1); y = p(:,2);
        rhs = x.*y;
    end
% gradient of exact solution
    function uprime = Du(p, ~)
        p = project(p);
        x = p(:,1); y = p(:,2); z = p(:,3);
        r = sqrt(x.^2+y.^2+z.^2);
        normal = [x./r, y./r, z./r];
        Du = zeros(size(x,1),3);
        Du(:,1) = y;
        Du(:,2) = x;
        uprime = Du - dot(Du, normal, 2).*normal;
        %uprime(:,4) = 0;
        %uprime = uprime;
    end
% Nueman boundary condition at t_0
    function rhs = g_N0(p)
        x = p(:,1); %y = p(:,2);
        rhs = zeros(size(x));
    end
% Nueman boundary condition at t_0
    function rhs = g_N1(p)
        x = p(:,1); %y = p(:,2);
        rhs = zeros(size(x));
    end
end
