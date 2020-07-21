clearvars();
close('all');
%clc();

%%
meshFile = 'meshes/commemi2d4.1';

mesh = mesh.getMesh('filename', meshFile, 'format', 'triangle', ...
    'shift', 0, 'scale', 1000, 'verbose', true);

%%
figure(1);
plot.plotMT('mesh', mesh, 'section', 'subdomains');

%%
freq = tools.pick(3, 1 / 1, 1 / 9, 1 / 100);
omega = 2 * pi * freq;
mu0 = pi * 4e-7;

sigma = 1.0 ./ [5, 1000, 10, 25, 1e9, 2.5, 2.5, 2.5, 2.5, 2.5, 1000];
mu = ones(size(sigma));
sigmaBCL = 1.0 ./ [25, 10, 1000, 5];
thkBCL = [0.5e3, 1.5e3, 23e3];
sigmaBCR = 1.0 ./ [25, 2.5, 1000, 5];
thkBCR = [0.5e3, 0.5e3, 24e3];

xobs = tools.asRow(-9e4:100:9e4);
obs = [xobs; zeros(size(xobs)) + 1];

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
