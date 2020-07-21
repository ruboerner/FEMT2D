clearvars();
close('all');
clc();

%%
meshFile = 'meshes/layered.1';

mesh = mesh.getMesh('filename', meshFile, ...
    'format', 'triangle', ...
    'shift', 0, 'scale', 1, 'verbose', true);

%%
figure(1);
plot.plotMT('mesh', mesh, 'section', 'subdomains');

%%
periods = logspace(0, 4, 41);
freq = 1 ./ periods;
omega = 2 * pi * freq;
mu0 = pi * 4e-7;

sigma = 1.0 ./ [100, 10, 100, 1e9];
mu = [1, 1, 1, 1];
sigmaBCL = 1.0 ./ [100, 10, 100];
muL = [1, 10, 1];
thkBCL = [18e3, 8e3];

xobs = tools.asRow(-1e4:1e3:1e4);
obs = [xobs; zeros(size(xobs)) + 0.1];

%%
fem = fe.FEMproblem('mesh', mesh, ...
    'elementtype', 'Lagrange', ...
    'order', 2, 'dimension', 2, ...
    'sigma', sigma, 'mu', mu, ...
    'application', 'MT', ...
    'polarization', 'both', ...
    'frequency', freq(1), ...
    'verbose', true);

%%
figure(2);
plot.plotMT('fem', fem, 'section', 'conductivity');

%%
fem = fe.getQ(fem, obs);
% fem = fe.getQfull(fem);

fem = fe.FEMassemble(fem, 'output', 'matrices', ...
    'verbose', false);

solution.rhoaxy = zeros(length(xobs), length(periods));
solution.rhoayx = zeros(length(xobs), length(periods));
solution.phixy = zeros(length(xobs), length(periods));
solution.phiyx = zeros(length(xobs), length(periods));


for k = 1:length(periods)
    freq = 1.0 / periods(k);
    
    fem = fe.removeDirichlet(fem, ...
        'sigmaBCL', sigmaBCL, ...
        'thicknessBCL', thkBCL, ...
        'frequency', freq);

    sol = fe.FEMsolve(fem, 'verbose', false);

    sol = mt.postProcessing(fem, sol);
    solution.rhoaxy(:, k) = sol.rhoaxy;
    solution.rhoayx(:, k) = sol.rhoayx;
    solution.phixy(:, k) = sol.phixy;
    solution.phiyx(:, k) = sol.phiyx;
end

%%
plot.plotMT('fem', fem, 'sol', solution, ...
    'sounding', 'rhoa+phase', ...
    'station', 11, 'periods', periods);
