function E = getE1dMT(f,rho,d,zi)
% E = getE1DMT(F,RHO,D,Z)
%
% Compute horizontal electric field E within and above a layered half-space.
% Plane wave excitation is assumed with magnetic field H = 1 A/m at z >= 0.
%
% Input:
%
% F:   Frequency
% RHO: NL vector of layer resistivites within conducting half-space
% D:   vector of layer thicknesses, length NL - 1
% Z:   NZ Vector of z coordinates where electric fields have to be computed.
%      Note that z axis is pointing downward.
%
% Output:
%
% E:   NZ vector of complex electric fields at depths Z.
% Note that E is normalized such that magnetic field H = 1 A/m at z = 0.
%
% Ralph-Uwe Boerner (2009)

% Homogeneous Halfspace in Comsol
% u = 0;
% rho = [0.01];
% d = [];

assert(nargin == 4, 'getE1DMT requires four input arguments.');

nl = length(rho);
rho = rho(:);
d = d(:);
% coordinates of interfaces
h = [0; cumsum(d)];
% Constants
sigair = 1e-14;
mu0 = 4e-7 * pi;
iwm = 1i * 2 * pi * f * mu0;
alpha = complex(zeros(nl, 1));
b = alpha;
aa = complex(zeros(nl - 1, 1));
nz = length(zi);
E = complex(zeros(nz, 1));

alpha = sqrt(iwm ./ rho);
if nl == 1
    c1 = iwm / alpha;
else
    alphad = alpha(1:(nl - 1)) .* d;
    talphad = tanh(alphad);
    % Compute admittance at surface of layer recursively
    b(nl) = alpha(nl);
    for nn = nl - 1:-1:1
        b(nn) = alpha(nn) * (b(nn + 1) + alpha(nn) * talphad(nn)) ./...
            (alpha(nn) + b(nn + 1) * talphad(nn));
    end
    % Impedance
    c1 = iwm ./ b(1);
    % Continuation from boundary to boundary
    for nn = 1:nl - 1
        aa(nn) = (b(nn) + alpha(nn)) / (b(nn + 1) + alpha(nn)) * ...
            exp( - alpha(nn) * d(nn));
    end
end

for ii = 1:nz
    z = zi(ii);
    if z >= 0
        if nl == 1
            a = exp(-alpha(nl) * z);
        else
            ind = find(z >= h, 1, 'last');
            if ind < nl
                a = prod(aa(1:ind - 1)) * 0.5 * (1 + b(ind) / alpha(ind)) * ...
                    (exp( - alpha(ind) * (z - h(ind))) - ...
                    (b(ind + 1) - alpha(ind)) / (b(ind + 1) + alpha(ind)) * ...
                    exp( - alpha(ind) * (d(ind) + h(ind + 1) - z)));
            else
                a = prod(aa) * exp( - alpha(ind) .* (z - h(ind)));
            end
        end
    else
        k0 = sqrt(iwm * sigair);
        pr = (c1 - iwm/k0) / (c1 + iwm/k0);
        ar = k0 * z;
        a = exp(-ar) * (1 + pr * exp(2 * ar)) * 1 ./ (1 + pr);
    end
    E(ii) = a * c1;
end
