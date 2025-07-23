function pde = TSElliott
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
pde = struct('phi',@phi,'gradient', @gradient, 'unitoutnormal',...
    @unitoutnormal, 'project', @project, 'initmesh',@initmesh,...
    'meancurvature', @meancurvature, 'exactu', @exactu, 'f', @f, 'Du',@Du, ...
    'g_N0',@g_N0, 'g_N1', @g_N1,'projectNT', @projectNT, 'DuDt', @DuDt);

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
        s = 1/4*x.^2 + y.^2 + 4*z.^2./(1+1/2*sin(pi*x)).^2-1;
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
        phix = x/2-4*pi*z.^2.*cos(pi*x)./(1+1/2*sin(pi*x)).^3;
        phiy = 2*y;
        phiz = 8*z./(1+1/2*sin(pi*x)).^2;
        n = [phix, phiy, phiz];
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

% project the surface using newtons method
    function node = projectNT(p)
        % projection function
        % Guess initial value
        node = p;
        valueAtNode = phi(node);
        normalAtNode = gradient(node);
        %vector = node - p;
        %e1 = normalAtNode./sqrt(dot(normalAtNode, normalAtNode,2));
        error=sqrt(valueAtNode.^2./(dot(normalAtNode,normalAtNode,2)));
        k=1;
        % Inital Guess
        bigX0 = [node, 2*valueAtNode./dot(normalAtNode,normalAtNode,2)];
        bigX = bigX0;
        x0 = bigX0(:,1); y0 = bigX0(:,2); z0 = bigX0(:,3);
        isNotOK = error > 1e-8;
        while  any(isNotOK) && k<200
            k=k+1;
            % Guess value
            x = bigX(isNotOK,1); y = bigX(isNotOK,2); z = bigX(isNotOK,3); lambda = bigX(isNotOK,4);
            % compute gradient for level set function
            phix = x/2-4*pi*z.^2.*cos(pi*x)./(1+1/2*sin(pi*x)).^3;
            phiy = 2*y;
            phiz = 8*z./(1+1/2*sin(pi*x)).^2;
            % compute hessian for level set function
            phixx = 1/2+4*pi^2*z.^2.*sin(pi*x)./(1+1/2*sin(pi*x)).^3 + ...
                6*pi.^2*z.^2.*cos(pi*x).^2./(1+1/2*sin(pi*x)).^4;
            phixy = zeros(size(y,1),1);
            phixz = -8*pi*z.*cos(pi*x)./(1+1/2*sin(pi*x)).^3;
            phiyy = 2*ones(size(y,1),1);
            phiyz = zeros(size(y,1),1);
            phizz = 8./(1+1/2*sin(pi*x)).^2;
            % assemble hessian for objective function F(X,lambda) =
            % |X-X0|^2 + lambda*phi(X) with X = [x,y,z]
            f11 = 2 + lambda.*phixx;
            f12 = lambda.*phixy;
            f13 = lambda.*phixz;
            f14 = phix;
            f21 = lambda.*phixy;
            f22 = 2 + lambda.*phiyy;
            f23 = lambda.*phiyz;
            f24 = phiy;
            f31 = lambda.*phixz;
            f32 = lambda.*phiyz;
            f33 = 2 + lambda.*phizz;
            f34 = phiz;
            f41 = phix;
            f42 = phiy;
            f43 = phiz;
            f44 = zeros(size(x));
            fhessian = [f11'; f21'; f31'; f41'; ...
                f12'; f22'; f32'; f42'; ...
                f13'; f23'; f33'; f43'; ...
                f14'; f24'; f34'; f44';];
            Mat = reshape(fhessian, 4, 4, size(x,1));
            % compute the gradient of objective function
            f1 = 2*(x-x0(isNotOK)) + lambda.*phix;
            f2 = 2*(y-y0(isNotOK)) + lambda.*phiy;
            f3 = 2*(z-z0(isNotOK)) + lambda.*phiz;
            f4 = phi(node(isNotOK,:));
            rhs = [ f1'; f2'; f3'; f4'];
            %rhs = reshape(fgrad, 4, 1, N);
            sol = MultiSolver(Mat, rhs);
            bigX(isNotOK,:) = bigX0(isNotOK,:) - sol';
            % compute error
            node(isNotOK,:) = bigX(isNotOK,[1,2,3]);
            valueAtNode = phi(node(isNotOK,:));
            normalAtNode = gradient(node(isNotOK,:));
            vector = node(isNotOK,:) - p(isNotOK);
            e1 = normalAtNode./sqrt(dot(normalAtNode, normalAtNode,2)) - ...
                vector./sqrt(dot(vector,vector,2));
            error=sqrt(valueAtNode.^2./(dot(normalAtNode, normalAtNode,2))+dot(e1,e1,2));
            isNotOK(isNotOK) = error > 1e-8;
        end
    end

% right hand side function
    function s = f(p, t)
        %p = project(p);
        x = p(:,1); y = p(:,2); z = p(:,3);
        phix = x/2-4*pi*z.^2.*cos(pi*x)./(1+1/2*sin(pi*x)).^3;
        phiy = 2*y;
        phiz = 8*z./(1+1/2*sin(pi*x)).^2;
        phixx = 1/2+4*pi^2*z.^2.*sin(pi*x)./(1+1/2*sin(pi*x)).^3 + ...
            6*pi.^2*z.^2.*cos(pi*x).^2./(1+1/2*sin(pi*x)).^4;
        phixy = zeros(size(y,1),1);
        phixz = -8*pi*z.*cos(pi*x)./(1+1/2*sin(pi*x)).^3;
        phiyy = 2*ones(size(y,1),1);
        phiyz = zeros(size(y,1),1);
        phizz = 8./(1+1/2*sin(pi*x)).^2;
        d = sqrt(phix.^2+phiy.^2+phiz.^2);
        hxx = phixx./d - phix.*(phix.*phixx+phiy.*phixy+phiz.*phixz)./d.^3;
        hyy = phiyy./d - phiy.*(phix.*phixy+phiy.*phiyy+phiz.*phiyz)./d.^3;
        hzz = phizz./d - phiz.*(phix.*phixz+phiy.*phiyz+phiz.*phizz)./d.^3;
        dn = hxx + hyy + hzz;
        s = 2*phix.*phiy./d.^2 + dn.*(y.*phix + x.*phiy)./d;
        s = s.*exp(t);
    end
% exact solution
    function rhs = exactu(p, t)
        %p = project(p);
        x = p(:,1); y = p(:,2);
        rhs = x.*y.*exp(t);
    end
% the gradient of exact solution in space
    function uprime = Du(p,t) % gradient of exact solution
        Du = zeros(size(p));
        Du(:,1)=p(:,2); Du(:,2)=p(:,1);Du(:,3)=0;
        normal = unitoutnormal(p);
        uprime = Du - dot(Du, normal, 2).*normal;
        %uprime(:,4) = 1;
        uprime = uprime.*exp(t);
    end
% the gradient of exact solution in time
function dudt = DuDt(p,t) % gradient of exact solution
        x = p(:,1); y = p(:,2);
        %uprime(:,4) = 1;
        dudt = x.*y.*exp(t);
    end
% Nueman boundary condition at t_0
    function rhs = g_N0(p)
        x = p(:,1); y = p(:,2);
        rhs = x.*y;
    end
% Nueman boundary condition at t_0
    function rhs = g_N1(p)
        x = p(:,1); y = p(:,2);
        rhs = x.*y*exp(1);
    end
end
