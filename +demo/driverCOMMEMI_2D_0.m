clearvars();
close('all');
clc();

%%
meshFile = 'meshes/commemi2d0.1';

mesh = mesh.getMesh('filename', meshFile, 'format', 'triangle', ...
    'shift', 0, 'scale', 1, 'verbose', true);

%%
figure(1);
plot.plotMT('mesh', mesh, 'section', 'subdomains');

%%
freq = 1 / 300;
omega = 2 * pi * freq;
mu0 = pi * 4e-7;

sigma = [1e1, 1e-1, 1e-9, 1, 0.5];
mu = ones(size(sigma));
sigmaBCL = [0.1, 1e1];
thkBCL = 5e4;
sigmaBCR = [0.5, 1e1];
thkBCR = 5e4;

xobs = tools.asRow(-4e4:250:4e4);
obs = [xobs; zeros(size(xobs)) + 1];

%%
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

fem = fe.removeDirichlet(fem, 'sigmaBCL', sigmaBCL, 'sigmaBCR', sigmaBCR, ...
    'thicknessBCL', thkBCL, 'thicknessBCR', thkBCR, 'frequency', freq);

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

%%
figure(6);
plot.plotMT('fem', fem, 'sol', sol, 'sounding', 'rhoa+phase', ...
    'station', 9);