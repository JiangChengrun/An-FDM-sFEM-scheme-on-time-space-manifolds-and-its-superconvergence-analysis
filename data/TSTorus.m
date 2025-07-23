function pde = TSTorus
%% TORUSPOISSONDATA Data for Laplace-beltrami operator
%   PDE = TORUSPOISSONDATA
%
%   Copyright (C) Hailong Guo
%   04/02/2015

%% Initalize
R = 4;
r = 1;

%% Main data structure
pde = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'initmesh',@initmesh,...
    'meancurvature', @meancurvature, 'exactu', @exactu, 'f', @f, 'Du',@Du, ...
    'g_N0',@g_N0, 'g_N1', @g_N1,'projectNT', @projectNT, 'DuDt', @DuDt);

%% Subfunctions
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
        s = x.^2 + y.^2 - 2*R*sqrt(x.^2 + y.^2) + R^2 + z.^2 - r^2;
    end

% unit out normal vector
    function un = unitoutnormal(p)
        n = gradient(p);
        nn = sqrt(sum(n.^2,2));
        un = n./repmat(nn,1,3);
    end

% gradient of level set function
    function n = gradient(p)
        % gradient of phi
        x = p(:,1); y = p(:,2); z = p(:,3);
        ux = 2*x - 2*R*x./sqrt(x.^2+y.^2);
        uy = 2*y - 2*R*y./sqrt(x.^2+y.^2);
        uz = 2*z;
        n = [ux, uy, uz];
    end

% project point on to the surface, this from is part of iFEM by Long Chen
    function [node,dist] = project(p)
        % projection function
        s = sign(phi(p));
        node = p;
        normalAtNode = gradient(node);
        valueAtNode = phi(node);
        node = node - valueAtNode*ones(1,3).*normalAtNode./(dot(normalAtNode,normalAtNode,2)*ones(1,3));
        vector = (-s)*ones(1,3).*(node - p);
        d = s.*sqrt(dot(vector,vector,2));
        normalAtNode = gradient(node);
        node = p - d*ones(1,3).*normalAtNode./(sqrt(dot(normalAtNode,normalAtNode,2))*ones(1,3));
        valueAtNode = phi(node);
        normalAtNode = gradient(node);
        vector = (-s)*ones(1,3).*(node - p);
        d = s.*sqrt(dot(vector,vector,2));
        e1 = normalAtNode./(sqrt(dot(normalAtNode, normalAtNode,2))*ones(1,3))-vector./(sqrt(dot(vector,vector,2))*ones(1,3));
        error=sqrt(valueAtNode.^2./(dot(normalAtNode, normalAtNode,2))+dot(e1,e1,2));
        k=1;
        while max(abs(error)) > 1e-6 && k<200
            k=k+1;
            node = node - valueAtNode*ones(1,3).*normalAtNode./(dot(normalAtNode,normalAtNode,2)*ones(1,3));
            vector = -s*ones(1,3).*(node - p);
            d = s.*sqrt(dot(vector,vector,2));
            normalAtNode = gradient(node);
            node = p - d*ones(1,3).*normalAtNode./(sqrt(dot(normalAtNode,normalAtNode,2))*ones(1,3));
            valueAtNode = phi(node);
            normalAtNode = gradient(node);
            vector = (-s)*ones(1,3).*(node - p);
            d = s.*sqrt(dot(vector,vector,2));
            e1 = normalAtNode./(sqrt(dot(normalAtNode, normalAtNode,2))*ones(1,3))-vector./(sqrt(dot(vector,vector,2))*ones(1,3));
            error=sqrt(valueAtNode.^2./(dot(normalAtNode, normalAtNode,2))+dot(e1,e1,2));
        end
        if nargout == 2
            dist = d;
        end
    end
    function rhs = f(p,t) % load data (right hand side function)
        x = p(:,1); y = p(:,2); z = p(:,3);
        rhs = (2*(x-y).*(sqrt(x.^2+y.^2)-4).*(-10*y.^2+x.^2.*(-10+sqrt(x.^2+y.^2))+ ...
            32*(-1+sqrt(x.^2+y.^2))-2*z.^2+sqrt(x.^2+y.^2).*(y.^2+z.^2)))./ ...
            ((x.^2 + y.^2).*(16 + x.^2 + y.^2 - 8*sqrt(x.^2 + y.^2) + z.^2).^2);
        rhs = rhs.*exp(t) ;
    end
% exact solution
    function rhs = exactu(p, t)
        p = project(p);
        x = p(:,1); y = p(:,2);
        rhs = (x-y).*exp(t);
    end
    function uprime = Du(p,t) % gradient of exact solution
        Du = zeros(size(p));
        Du(:,1)= 1; Du(:,2)=-1;Du(:,3)=0;
        normal = unitoutnormal(p);
        uprime = Du - dot(Du, normal, 2).*normal;
        %uprime(:,4) = 1;
        uprime = uprime.*exp(t);
    end
    function dudt = DuDt(p,t) % gradient of exact solution
        x = p(:,1); y = p(:,2);
        %uprime(:,4) = 1;
        dudt = (x-y).*exp(t);
    end
% Nueman boundary condition at t_0
    function rhs = g_N0(p)
        x = p(:,1); y = p(:,2);
        rhs = x-y;
    end
% Nueman boundary condition at t_0
    function rhs = g_N1(p)
        x = p(:,1); y = p(:,2);
        rhs = (x-y)*exp(1);
    end
end
