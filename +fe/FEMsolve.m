function sol = FEMsolve(fem, varargin)
%sol = FEMsolve(fem)

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