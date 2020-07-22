clearvars();
close('all');

%%
meshFile = 'meshes/model_with_topography.1';

mesh = mesh.getMesh('filename', meshFile, 'format', 'triangle', ...
    'shift', 0, 'scale', 1, 'verbose', true);

%%
figure(1);
plot.plotMT('mesh', mesh, 'section', 'subdomains');

%%
freq = 1 ./ 10;
omega = 2 * pi * freq;
mu0 = pi * 4e-7;

sigma = [1e-3, 1e-3, 1e-9, 1e-3, 1e-3, 1e-3];
mu = ones(size(sigma));
sigmaBC = 1e-3;
thkBC = [];

xobs = tools.asRow(-3.9e4:250:3.9e4);
obs = [xobs; -exp(-1e-7*xobs.^2) + 1];

fem = fe.FEMproblem('mesh', mesh, ...
    'elementtype', 'Lagrange', ...
    'order', 2, 'dimension', 2, ...
    'sigma', sigma, 'mu', mu, ...
    'application', 'MT', ...
    'polarization', 'both', ...
    'frequency', freq, ...
    'verbose', true);

%%
figure(2);
plot.plotMT('fem', fem, 'section', 'conductivity');

%%
fem = fe.getQ(fem, obs);
fem = fe.getQfull(fem);

fem = fe.FEMassemble(fem, 'output', 'matrices', 'verbose', true);

fem = fe.removeDirichlet(fem, 'sigmaBCL', sigmaBC, 'sigmaBCR', sigmaBC, ...
    'thicknessBCL', thkBC, 'thicknessBCR', thkBC, 'frequency', freq);

sol = fe.FEMsolve(fem, 'verbose', true);

sol = mt.postProcessing(fem, sol);

%%
figure(3);
plot.plotMT('fem', fem, 'sol', sol, 'obs', obs, 'profile', 'rhoa+phase');

%%
figure(4);
plot.plotMT('fem', fem, 'sol', sol, 'section', 'Jx');

%%
figure(5);
plot.plotMT('fem', fem, 'sol', sol, 'obs', obs, 'profile', 'tipper');
