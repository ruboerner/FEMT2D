function [Phi, GradPhi] = getBasisLagrange(x, y, order, dim)
%getBasisLagrange Compute finite element basis and its gradient for first
%and second order Lagrange (H1) finite elements at point (x,y)
%
% [Phi, gradPhi] = getBasisLagrange(x, y, order, dim)
%
% getBasisLagrange computes the Lagrange basis and its gradient
% within a triangular reference element connecting the points
% (0,0), (1,0), and (0,1).
%
% (C) 2020 Ralph-Uwe Börner
%

if nargin < 4
    dim = 2;
end

if nargin < 3
    order = 1;
end

switch dim
    case 2
        GradLambda = [...
            1 0 -1; ...
            0 1 -1];
        switch order
            case 1
                Phi = [x, y, 1 - x - y];
                GradPhi = GradLambda;
            case 2
                Lambda = [x, y,  1 - x - y];
                
                Phi(1) = Lambda(1)*(2*Lambda(1)-1);
                Phi(2) = Lambda(2)*(2*Lambda(2)-1);
                Phi(3) = Lambda(3)*(2*Lambda(3)-1);
                Phi(4) = 4*Lambda(2)*Lambda(3);
                Phi(5) = 4*Lambda(1)*Lambda(3);
                Phi(6) = 4*Lambda(1)*Lambda(2);
                
                GradPhi(:,1) = (4*Lambda(1)-1) * GradLambda(:,1);
                GradPhi(:,2) = (4*Lambda(2)-1) * GradLambda(:,2);
                GradPhi(:,3) = (4*Lambda(3)-1) * GradLambda(:,3);
                GradPhi(:,4) = 4*(Lambda(2)*GradLambda(:,3)+Lambda(3)*GradLambda(:,2));
                GradPhi(:,5) = 4*(Lambda(1)*GradLambda(:,3)+Lambda(3)*GradLambda(:,1));
                GradPhi(:,6) = 4*(Lambda(1)*GradLambda(:,2)+Lambda(2)*GradLambda(:,1));
            otherwise
                error('Higher order elements not implemented.');
        end
        
    otherwise
        error('3-D not implemented.');
end
