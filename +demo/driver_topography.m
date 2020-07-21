clearvars();
close('all');

%%
meshFile = 'meshes/model_with_topography.1';

mesh = mesh.getMesh('filename', meshFile, 'format', 'triangle', ...
    'shift', 0, 'scale', 1, 'verbose', true);

%%
plot.plotSubdomains(mesh, 'meshcolor', 'none');
xlabel('y in m');
ylabel('z in m');
axis equal tight ij

%%
freq = 1;
omega = 2 * pi * freq;
mu0 = pi * 4e-7;

sigma = [1e-9, 2e-3, 2e-3, 5e-0, 2e-3, 6e-0];
sigma = [1e-3, 1e-3, 1e-3, 1e-9, 1e-3, 1e-1];
sigma = [1e-3, 1e-3, 1e-9, 1e-3, 1e-3, 1e-3];
mu = ones(size(sigma));
mu(5) = 1.0;
sigmaBC = 1e-3;
thkBC = [];

% xobs = asRow(0:1:30);
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
figure(1);
plot.patchplotConst(fem.mesh.node, ...
    [fem.mesh.tri2node; tools.asRow(fem.sigma)], @log10, [0.9 0.9 0.9]);
xlabel('y in m');
ylabel('z in m');
colormap(gca, jet(256));
h = colorbar();
axis equal tight ij
ylabel(h, 'log_{10}(\sigma/({S/m}))');
% 
% title('Conductivity');

%%
fem = fe.getQ(fem, obs);
fem = fe.getQfull(fem);

fem = fe.FEMassemble(fem, 'output', 'matrices', 'verbose', true);

fem = fe.removeDirichlet(fem, 'sigmaBCL', sigmaBC, 'sigmaBCR', sigmaBC, ...
    'thicknessBCL', thkBC, 'thicknessBCR', thkBC, 'frequency', freq);

sol = fe.FEMsolve(fem, 'verbose', true);

iw = -2i * pi * freq;
Ex = fem.Q.Qe   * sol.ue;
Ey = fem.Q.QhEy * sol.uh;
Ez = fem.Q.QhEz * sol.uh;
Hx = fem.Q.Qh   * sol.uh;
Hy = fem.Q.QeHy * sol.ue / iw;
Hz = fem.Q.QeHz * sol.ue / iw;

T = Hz ./ Hy;

Zxy = Ex ./ Hy;
rhoaxy = abs(Zxy).^2 / mu0 / omega;
phixy = atan(imag(Zxy) ./ real(Zxy)) * 180 / pi;

Zyx = Ey ./ Hx;
rhoayx = abs(Zyx).^2 / mu0 / omega;
phiyx = atan(imag(Zyx) ./ real(Zyx)) * 180 / pi;

%%
figure(2);
subplot(2, 1, 1);
semilogy(xobs, rhoaxy, '.-', xobs, rhoayx, '.-');
ylim([100 10000]);
grid();
xlabel('y in km');
legend('E-pol', 'H-Pol');
subplot(2, 1, 2);
plot(xobs, phixy, '.-', xobs, phiyx, '.-');
ylim([0 90]);
grid();
xlabel('y in km');
legend('E-pol', 'H-Pol');

%%
figure(3);
plot.patchplotConst(fem.mesh.node, ...
    [fem.elem2dofs(1:3, :); tools.asRow(abs(fem.Q.QJefull * sol.ue))], ...
    @(x) x, 'none');
colormap(gca, jet(256));
h = colorbar();
ylabel(h, '|J| in A/m^2');
axis equal tight ij

%%
figure(4);
plot(xobs, real(T), 'r--', xobs, imag(T), 'b--');
grid();
xlabel('y in km');
ylabel('magn. transfer function');
legend('Re', 'Im');