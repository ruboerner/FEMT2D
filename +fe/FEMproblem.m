function fem = FEMproblem(varargin)
%FEMproblem provides struct fem
%
% fem = FEMproblem() computes all FE related quantities necessary for the
% solution of the 2D MT problem.
%
% Input:
% ======
%
% FEMproblem() accepts the following keyword-value
% pairs:
%
% | 'mesh'         | mesh struct from call to getMesh()
% | 'order'        | FE order 1 or 2 ((linear or quadratic)
% | 'polarization' | 'epol', 'hpol', or 'both'
% | 'sigma'        | array of conductivities given in the order of
% |                | subdomains
% | 'mu'           | likewise for *relative* magnetic permeability
% | 'frequency'    | frequency in Hz, scalar
% | 'verbose'      | if true, FEMproblem() provides more information during execution
%
% Output:
% =======
%
% FEMproblem() returns the MATLAB struct fem that holds all information necessary for later assembly
% of FE matrices etc.
%
% (C) 2020 Ralph-Uwe Börner


p = inputParser;

addParameter(p, 'mesh', struct([]), @isstruct);
addParameter(p, 'order', 1, @isnumeric);
addParameter(p, 'dimension', 2, @isnumeric);
addParameter(p, 'application', 'none', @ischar);
addParameter(p, 'polarization', 'none', @ischar);
addParameter(p, 'elementtype', 'Lagrange', @ischar);
addParameter(p, 'sigma', [], @isnumeric);
addParameter(p, 'mu', [], @isnumeric);
addParameter(p, 'frequency', [], @isnumeric);
addParameter(p, 'time', [], @isnumeric);
addParameter(p, 'wavenumber', [], @isnumeric);
addParameter(p, 'verbose', true, @islogical);

parse(p, varargin{:});

verbose = p.Results.verbose;
fem.mesh = p.Results.mesh;
fem.order = p.Results.order;
fem.dimension = p.Results.dimension;
fem.application = p.Results.application;

fem.app.frequency = p.Results.frequency;
fem.app.time = p.Results.time;
fem.app.wavenumber = p.Results.wavenumber;


if strcmp(fem.application, 'MT')
    assert(strcmp(p.Results.polarization, 'epol') ...
        | strcmp(p.Results.polarization, 'hpol') ...
        | strcmp(p.Results.polarization, 'both'));
end
fem.polarization = p.Results.polarization;


fem.elementtype = p.Results.elementtype;

sigma = p.Results.sigma;
assert(length(sigma) == length(unique(fem.mesh.tri2subdomain)), ...
    'FEMproblem::Number of entries in sigma must match number of subdomains.');

mu = p.Results.mu;
if isempty(mu)
    mu = ones(size(sigma));
else
    assert(length(mu) == length(unique(fem.mesh.tri2subdomain)), ...
        'FEMproblem::Number of entries in mu must match number of subdomains.');
end

sigma_ = zeros(fem.mesh.nt, 1);
mu_ = zeros(fem.mesh.nt, 1);

for k = 1:fem.mesh.nt
    sigma_(k) = sigma(fem.mesh.tri2subdomain(k));
    mu_(k) = 4e-7 * pi * mu(fem.mesh.tri2subdomain(k));
end

fem.sigma = sigma_(:);
fem.mu = mu_(:);

fem.stripair = false;
fem.istrip = 0;

if strcmp(fem.polarization, 'hpol') || strcmp(fem.polarization, 'both')
    fem.stripair = true;
    %     fem.istrip = find(fem.sigma == min(sigma));
    fem.istrip = find(sigma==min(sigma));
end

fem.coorddofs = zeros(2, fem.mesh.np + (fem.order - 1) * fem.mesh.ne);

switch fem.order
    case 1
        fem.elem2dofs = fem.mesh.tri2node;
        fem.ndofs = fem.mesh.np;
        fem.bddofs = fem.mesh.bdNodes;
        C = unique(reshape(fem.elem2dofs(1:3, :), [], 1));
        fem.coorddofs = fem.mesh.node(:, C);
    case 2
        fem.elem2dofs = [fem.mesh.tri2node; fem.mesh.np + fem.mesh.tri2edge];
        fem.ndofs = fem.mesh.np + fem.mesh.ne;
        fem.bddofs = [fem.mesh.bdNodes; fem.mesh.np + fem.mesh.bdEdges];
        C = unique(reshape(fem.elem2dofs(1:3, :), [], 1));
        fem.coorddofs = fem.mesh.node(:, C);
        for k = 1:fem.mesh.ne
            fem.coorddofs(1, fem.mesh.np + k) = ...
                mean(fem.mesh.node(1, fem.mesh.edge2node(:, k)));
            fem.coorddofs(2, fem.mesh.np + k) = ...
                mean(fem.mesh.node(2, fem.mesh.edge2node(:, k)));
        end
end
 
xleftbd = min(fem.coorddofs(1,fem.bddofs));
xrightbd = max(fem.coorddofs(1,fem.bddofs));

ztopbd = min(fem.coorddofs(2,fem.bddofs));
zbotbd = max(fem.coorddofs(2,fem.bddofs));

idleft = find(abs(fem.coorddofs(1, fem.bddofs) - xleftbd) < sqrt(eps));
idright = find(abs(fem.coorddofs(1, fem.bddofs) - xrightbd) < sqrt(eps));

idtop = find(abs(fem.coorddofs(2, fem.bddofs) - ztopbd) < sqrt(eps));
idbot = find(abs(fem.coorddofs(2, fem.bddofs) - zbotbd) < sqrt(eps));


fem.bddofsleft  = fem.bddofs(idleft);    %#ok
fem.bddofsright = fem.bddofs(idright);   %#ok
fem.bddofstop   = fem.bddofs(idtop);     %#ok
fem.bddofsbot   = fem.bddofs(idbot);     %#ok

indexDirichlet = fem.bddofs;
indexInner = setdiff(1:fem.ndofs, indexDirichlet);
ndofDirichlet = length(indexDirichlet);
ndofInner = length(indexInner);

fem.Compl = sparse(...
    indexDirichlet, ...
    1:ndofDirichlet, ...
    ones(ndofDirichlet, 1), ...
    fem.ndofs, ndofDirichlet);
fem.Null = sparse(...
    indexInner, ...
    1:ndofInner, ...
    ones(ndofInner, 1), ...
    fem.ndofs, ndofInner);

if fem.stripair
    indexDOFS = unique([indexDirichlet; fe.finddofsinsubdomain(fem, fem.istrip)]);
    fem.stripdofs = indexDOFS;
    indexOther = setdiff(1:fem.ndofs, indexDOFS);
    ndofDOFS = length(indexDOFS);
    ndofOther = length(indexOther);
    fem.ComplStrip = sparse(...
        indexDOFS, ...
        1:ndofDOFS, ...
        ones(ndofDOFS, 1), ...
        fem.ndofs, ndofDOFS);
    fem.NullStrip = sparse(...
        indexOther, ...
        1:ndofOther, ...
        ones(ndofOther, 1), ...
        fem.ndofs, ndofOther);
end

if verbose
    fprintf('\nFE summary:\n');
    fprintf('===========\n');
    fprintf('  Application mode: %s\n', fem.application);
    if strcmp(fem.application, 'MT')
        fprintf('  MT polarization(s): %s\n', fem.polarization);
    end
    fprintf('  Element type: %s\n', fem.elementtype);
    fprintf('  Dimension: %d\n', fem.dimension);
    fprintf('  FE order: %d\n', fem.order);
    fprintf('  DOFS: %d\n', fem.ndofs);
end