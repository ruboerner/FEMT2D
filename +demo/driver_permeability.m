clearvars();
close('all');
%clc();

%%
meshFile = 'meshes/permeability.1';

mesh = mesh.getMesh('filename', meshFile, 'format', 'triangle', ...
    'shift', 0, 'scale', 1, 'verbose', true);

%%
figure(1);
plot.plotMT('mesh', mesh, 'section', 'subdomains');

%%
freqs = 2.^(-13:3);

freq = freqs(1);

omega = 2 * pi * freq;
mu0 = pi * 4e-7;

sigma = 1.0 ./ [100, 500, 4000, 25, 1e9, 10, 25, 25, 25, 25, 4000];
mu = ones(size(sigma));
mu([8, 10]) = 1;
sigmaBCL = 1.0 ./ [25, 4000, 500, 100];
thkBCL = [1e3, 99e3, 100e3];
sigmaBCR = sigmaBCL;
thkBCR = thkBCL;

xobs = tools.asRow(-180e3:1000:320e3);
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
