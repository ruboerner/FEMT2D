function fem = getQfull(fem)
%fem = getQfull(fem)
%

nD = (fem.dimension + 1) + (fem.order - 1) * fem.dimension * (fem.dimension + 1) / 2;

t = fem.elem2dofs;
p = fem.mesh.node;
ndofs = fem.ndofs;
nt = fem.mesh.nt;

Qe = zeros(nD, nt);
QJe = zeros(nD, nt);
QeHy = zeros(nD, nt);
QeHz = zeros(nD, nt);
Qh = zeros(nD, nt);
QhEy = zeros(nD, nt);
QhEz = zeros(nD, nt);


[phi_hat, grad] = fe.getBasisLagrange(1/3, 1/3, fem.order);

for indEl = 1:nt
    B_k = fem.mesh.Bk(:, :, indEl);
    b_k = p(1:2, t(3, indEl));
    r = B_k * [1/3; 1/3] + b_k;
    B = B_k' \ grad;
    
    if strcmp(fem.polarization, 'epol') || strcmp(fem.polarization, 'both')
        Qe(:, indEl) = phi_hat';
        QJe(:, indEl) = fem.sigma(indEl) * Qe(:, indEl);
        QeHy(:, indEl) = -B(2, :) ./ fem.mu(indEl);
        QeHz(:, indEl) =  B(1, :) ./ fem.mu(indEl);
    end
    
    if strcmp(fem.polarization, 'hpol') || strcmp(fem.polarization, 'both')
        Qh(:, indEl) = phi_hat';
        QhEy(:, indEl) =  B(2, :) ./ fem.sigma(indEl);
        QhEz(:, indEl) = -B(1, :) ./ fem.sigma(indEl);
    end
    
    if strcmp(fem.application, 'axisymm')
        Qe(:, indEl) = phi_hat';
        QJe(:, indEl) = fem.sigma(indEl) * Qe(:, indEl);
        QeHy(:, indEl) = -B(2, :);
        QeHz(:, indEl) = B(1, :) + phi_hat / r(1);
    end
end

Ci = repmat(tools.asRow(1:nt), nD, 1);

if strcmp(fem.polarization, 'epol') || strcmp(fem.polarization, 'both')
    fem.Q.Qefull = sparse(Ci(:), t(:), Qe(:), nt, ndofs);
    fem.Q.QJefull = sparse(Ci(:), t(:), QJe(:), nt, ndofs);
    fem.Q.QeHyfull = sparse(Ci(:), t(:), QeHy(:), nt, ndofs);
    fem.Q.QeHzfull = sparse(Ci(:), t(:), QeHz(:), nt, ndofs);
end
if strcmp(fem.polarization, 'hpol') || strcmp(fem.polarization, 'both')
    fem.Q.Qhfull = sparse(Ci(:), t(:), Qh(:), nt, ndofs);
    fem.Q.QhEyfull = sparse(Ci(:), t(:), QhEy(:), nt, ndofs);
    fem.Q.QhEzfull = sparse(Ci(:), t(:), QhEz(:), nt, ndofs);
end
if strcmp(fem.application, 'TEM')
    fem.Q.Qefull = sparse(Ci(:), t(:), Qe(:), nt, ndofs);
    fem.Q.QJefull = sparse(Ci(:), t(:), QJe(:), nt, ndofs);
    fem.Q.QeHyfull = sparse(Ci(:), t(:), QeHy(:), nt, ndofs);
    fem.Q.QeHzfull = sparse(Ci(:), t(:), QeHz(:), nt, ndofs);
end
