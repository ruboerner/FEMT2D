function sol = FEMsolve(fem, varargin)
%FEMsolve solve system of linear equations for the MT problem
%
% sol = FEMsolve(fem) computes the numerical solution of the
% linear systems arising in the 2D MT problem.
%
% sol = FEMsolve(fem, 'verbose', true) computes the solution and displays
% additional timing information.
%
% Input:
% ======
%
% fem struct resulting from a call to FEMproblem()
%
% Output:
% =======
% The struct sol contains the numerical solution (i.e., the degrees of
% freedom) out of which the MT fields can be computed in the
% post-processing stage. More precisely, the struct sol has the fields
%
% sol.ue for E-polarisation
% sol.uh for H-polarisation
%
% (C) 2020 Ralph-Uwe Börner

p = inputParser;
addParameter(p, 'verbose', true, @islogical);
% addParameter(p, 'frequency', [], @isnumerical);

parse(p, varargin{:});

verbose = p.Results.verbose;

if verbose
    fprintf('\nFE solver:');
    fprintf('\n==========\n');
    fprintf('  Linear system solve...');
    t0 = tic();
end

if ~strcmp(fem.application, 'MT')
    sol.ue = fem.uDir + fem.Null * (fem.Ac \ fem.fc);
end

if strcmp(fem.polarization, 'epol') || strcmp(fem.polarization, 'both')
    sol.ue = fem.uDir + fem.Null * (fem.Ac \ fem.fc);
end
if strcmp(fem.polarization, 'hpol') || strcmp(fem.polarization, 'both')
    sol.uh = fem.uDirh + fem.NullStrip * (fem.Ach \ fem.fch);
end

if verbose
    tSolve = toc(t0);
    fprintf('done (%.2f s).\n', tSolve);
end