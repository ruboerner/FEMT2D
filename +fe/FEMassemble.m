function fem = FEMassemble(fem, varargin)

p = inputParser;
addParameter(p, 'output', 'matrices', @ischar);
addParameter(p, 'verbose', true, @islogical);


parse(p, varargin{:});

outputformat = p.Results.output;
verbose = p.Results.verbose;

% Size of element matrices is nD x nD:
%
nD = (fem.dimension + 1) + (fem.order - 1) * fem.dimension * (fem.dimension + 1) / 2;

if verbose
    fprintf('\nAssembly:\n');
    fprintf('=========\n');
    fprintf('  Size of element matrices: %d x %d\n', nD, nD);
end

% Quadrature nodes and weights
switch fem.order
    case 1
        xi    = 1/3;
        eta   = 1/3;
        w     = 1/2;
        M_ref = 1 / 24 * ...
            [2 1 1; ...
            1 2 1; ...
            1 1 2];
    case 2
        xi    = [1/2, 1/2,   0];
        eta   = [  0, 1/2, 1/2];
        w     = [1/6, 1/6, 1/6];
        M_ref = [
            1/60 -1/360 -1/360     0     0 -1/90;
            -1/360   1/60 -1/360     0 -1/90     0;
            -1/360 -1/360   1/60 -1/90     0     0;
            0      0  -1/90  4/45  2/45  2/45;
            0  -1/90      0  2/45  4/45  2/45;
            -1/90      0      0  2/45  2/45  4/45
            ];
end

nQuad = length(w);

% Evaluate basis function and gradient at quadrature nodes xi and eta
%
phi = zeros(nD, nQuad);
gradphi = zeros(2, nD, nQuad);

% Evaluate basis function and its gradient at quadrature nodes within
% reference element
%
for ii = 1:length(w)
    [phi(:, ii), gradphi(:, :, ii)] = fe.getBasisLagrange(xi(ii), eta(ii), fem.order);
end

% Number of triangular elements
%
nEl = fem.mesh.nt;

% Allocate storage for all element matrices in nD x nD x nEl arrays
%
nDOFS = fem.ndofs;
K = zeros(nD, nD, nEl);
M = zeros(nD, nD, nEl);
if strcmp(fem.application, 'csem')
    L = zeros(nD, nD, NEl);
end


Ci = repmat(reshape(fem.elem2dofs, nD, 1, []), [1, nD, 1]);
Cj = repmat(reshape(fem.elem2dofs, 1, nD, []), [nD, 1, 1]);

% Assemble element matrices
%
if verbose
    fprintf('  Assembly of element matrices...');
    t0 = tic();
end

for i = 1:nEl
    K_k = zeros(nD, nD);
    F_k = zeros(nD, nD);
    if strcmp(fem.application, 'csem')
        L_k = zeros(nD, nD);
    end
    if strcmp(fem.application, 'axisymm')
        b_k = fem.mesh.node(:, fem.elem2dofs(3, i));
    end
    d = abs(fem.mesh.detBk(i));
    Binv = fem.mesh.Binv(:, :, i);
    B_k = fem.mesh.Bk(:, :, i);
    for k = 1:nQuad
        B = Binv * gradphi(:, :, k);
        if strcmp(fem.application, 'axisymm')
            pt = B_k * [xi(k); eta(k)] + b_k;
            K_k = K_k + w(k) * (B' * B) * abs(pt(1));
        elseif strcmp(fem.application, 'csem')
            L_k = L_k + w(k) * (B' * [0 -1; 1 0] * B);
        else
            K_k = K_k + w(k) * (B' * B);
        end
    end
    
    M_k = M_ref * d;
    K_k = K_k   * d;
    if strcmp(fem.application, 'axisymm')
        rmid = abs(B_k * [1/3; 1/3] + b_k);
        M_k = M_ref * d * rmid(1);
        F_k = M_ref * d / rmid(1);
    end
    M(:, :, i) = M_k;
    K(:, :, i) = K_k + F_k;
    if strcmp(fem.application, 'csem')
        L(:, :, i) = L_k * d;
    end
end

if verbose
    tAssembly = toc(t0);
    fprintf('done (%.2f s).\n', tAssembly);
end

if strcmp(outputformat, 'matrices')
    if verbose
        fprintf('  Assembly of system matrices...');
        t0 = tic();
    end
    
    % Multiply with parameter
    %
    fem.Ke = sparse(Ci(:), Cj(:), tools.asColumn(bsxfun(@times, K, reshape(1 ./ fem.mu, 1, 1, []))), nDOFS, nDOFS);
    fem.Me = sparse(Ci(:), Cj(:), tools.asColumn(bsxfun(@times, M, reshape(fem.sigma, 1, 1, []))), nDOFS, nDOFS);
    
    fem.Kh = sparse(Ci(:), Cj(:), tools.asColumn(bsxfun(@times, K, reshape(1 ./ fem.sigma, 1, 1, []))), nDOFS, nDOFS);
    fem.Mh = sparse(Ci(:), Cj(:), tools.asColumn(bsxfun(@times, M, reshape(fem.mu, 1, 1, []))), nDOFS, nDOFS);
end

if verbose
    tAssembly = toc(t0);
    fprintf('done (%.2f s).\n', tAssembly);
end

