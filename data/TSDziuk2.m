function pde = TSDziuk2
%% TORUSPOISSONDATA Data for Laplace-beltrami operator
%   PDE = TORUSPOISSONDATA
%
%   Copyright (C) Hailong Guo
%   04/02/2015

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
        s =  (p(:,1) - p(:,3) .^ 2) .^ 2 + p(:,2) .^ 2 + p(:,3) .^ 2 - 1;
    end

% unit out normal vector
    function un = unitoutnormal(p)
        n = gradient(p);
        nn = sqrt(sum(n.^2,2));
        un = n./repmat(nn,1,3);
        %         s = 1 + 4 .* p(:,3) .^ 6 + (-8 .* p(:,1) + 4) .* p(:,3) .^ 4 + (-4 .* p(:,1) + 4 .* p(:,1) .^ 2) .* p(:,3) .^ 2;
        %         n = [(s .^ (-0.1e1 ./ 0.2e1) .* (2 .* p(:,1) - 2 .* p(:,3) .^ 2)) ./ 0.2e1, s .^ (-0.1e1 ./ 0.2e1) .* p(:,2),...
        %             (s .^ (-0.1e1 ./ 0.2e1) .* (-4 .* (p(:,1) - p(:,3) .^ 2) .* p(:,3) + 2 .* p(:,3))) ./ 0.2e1];
    end

% gradient of level set function
    function n = gradient(p)
        % gradient of phi
        n = [2 .* p(:,1) - 2 .* p(:,3) .^ 2, 2 .* p(:,2), ...
            -4 .* (p(:,1) - p(:,3) .^ 2) .* p(:,3) + 2 .* p(:,3)];
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
    function rhs = f(p,~) % load data (right hand side function)
        x = p(:,1); y = p(:,2); z = p(:,3);
        rhs = 2*y.*(2*x-z.^2).*(-x.^3 + z.^2 + 3*z.^4 + 3*z.^6 + ...
            y.^2.*(1 + 3*z.^2) + x.^2.*(1 + 5*z.^2) - ...
            x.*(y.^2 + 4.*z.^2 + 7*z.^4))./...
            (y.^2 + z.^2 + 5.*z.^4 + 4.*z.^6 + x.^2.*(1 + 4*z.^2) - ...
            2*x.*z.^2.*(3 + 4*z.^2)).^2 + ...
            8*y.*(x - z.^2)./(4*y.^2 + 4*(x - z.^2).^2 + ...
            (2*z - 4*z.*(x - z.^2)).^2);
        rhs = rhs + x.*y ;
    end
% exact solution
    function rhs = exactu(p, ~)
        p = project(p);
        x = p(:,1); y = p(:,2);
        rhs = x.*y;
    end
    function uprime = Du(p,~) % gradient of exact solution
        Du = zeros(size(p));
        Du(:,1)=p(:,2); Du(:,2)=p(:,1);Du(:,3)=0;
        normal = unitoutnormal(p);
        uprime = Du - dot(Du, normal, 2).*normal;
        %uprime(:,4) = 0;
        %uprime = uprime;
    end
function dudt = DuDt(p,t) % gradient of exact solution
        x = p(:,1); y = p(:,2);
        %uprime(:,4) = 1;
        dudt = 0*x.*y.*exp(t);
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