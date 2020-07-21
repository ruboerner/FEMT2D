function fem = getQ(fem, obs)
%function fem = getQ(fem, obs)

% dim = fem.dimension;
% order = fem.order;

% Identify elements to which points (x, y) belong ans store their indices
% in array indEl
%
indEl = tsearchn(fem.mesh.node(1:2, :)', fem.mesh.tri2node(1:3, :)', ...
    [tools.asColumn(obs(1, :)), tools.asColumn(obs(2, :))]);

nEl = length(indEl);

nD = (fem.dimension + 1) + (fem.order - 1) * fem.dimension * (fem.dimension + 1) / 2;

Qe = zeros(nD, nEl);
QeHy = zeros(nD, nEl);
QeHz = zeros(nD, nEl);
if strcmp(fem.polarization, 'hpol') || strcmp(fem.polarization, 'both')
    Qh = zeros(nD, nEl);
    QhEy = zeros(nD, nEl);
    QhEz = zeros(nD, nEl);
end

Ci = repmat(tools.asRow(1:nEl), nD, 1);
Cj = zeros(nD, nEl);

for k = 1:nEl
    B_k = fem.mesh.Bk(:, :, indEl(k));
    b_k = fem.mesh.node([1,2], fem.mesh.tri2node(3, indEl(k)));
    
    xi_hat = B_k \ ([obs(1, k); obs(2, k)] - b_k);
    [phi_hat, grad] = fe.getBasisLagrange(xi_hat(1), xi_hat(2), fem.order);
    
    B = B_k' \ grad;
    Cj(:, k) = tools.asColumn(fem.elem2dofs(:, indEl(k)));
    if strcmp(fem.polarization, 'epol') || strcmp(fem.polarization, 'both')
        Qe(:, k) = phi_hat;
        QeHy(:, k) = +B(2, :) ./ fem.mu(indEl(k));
        QeHz(:, k) = -B(1, :) ./ fem.mu(indEl(k));
    end
    if strcmp(fem.polarization, 'hpol') || strcmp(fem.polarization, 'both')
        Qh(:, k) = phi_hat;
        QhEy(:, k) = +B(2, :) ./ fem.sigma(indEl(k));
        QhEz(:, k) = -B(1, :) ./ fem.sigma(indEl(k));
    end
    if strcmp(fem.application, 'axisymm')
        Qe(:, k) = phi_hat;
        QeHy(:, k) = +B(2, :);
        QeHz(:, k) = -B(1, :) + phi_hat / obs(1, k);
    end
end

fem.Q.Qe = sparse(Ci(:), Cj(:), Qe(:), nEl, fem.ndofs);
fem.Q.QeHy = sparse(Ci(:), Cj(:), QeHy(:), nEl, fem.ndofs);
fem.Q.QeHz = sparse(Ci(:), Cj(:), QeHz(:), nEl, fem.ndofs);

if strcmp(fem.polarization, 'hpol') || strcmp(fem.polarization, 'both')
    fem.Q.Qh = sparse(Ci(:), Cj(:), Qh(:), nEl, fem.ndofs);
    fem.Q.QhEy = sparse(Ci(:), Cj(:), QhEy(:), nEl, fem.ndofs);
    fem.Q.QhEz = sparse(Ci(:), Cj(:), QhEz(:), nEl, fem.ndofs);
end

end