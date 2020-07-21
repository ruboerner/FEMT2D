clearvars();
close('all');
%clc();

%%
meshFile = 'meshes/weaver.1';

mesh = mesh.getMesh('filename', meshFile, 'format', 'triangle', ...
    'shift', 0, 'scale', 1, 'verbose', true);

%%
plot.plotSubdomains(mesh, 'meshcolor', 'none');
xlabel('y in m');
ylabel('z in m');
axis equal tight ij

%%
periods = 100.0;
freq = 1 ./ periods;
omega = 2 * pi * freq;
mu0 = pi * 4e-7;

sigma = [1e3, 1e-1, 1e-9, 1e-2];
mu = ones(size(sigma));
sigmaBCL = [0.1, 1e3];
thkBCL = 200e3;
sigmaBCR = [0.01, 1e3];
thkBCR = 200e3;

xobs = tools.asRow([-5e4 -1e4 -5e3:1e3:-1e3 -1 +1 1e3:1e3:5e3 1e4 5e4]);
obs = [xobs; zeros(size(xobs)) + 0.00001];

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

%%
fem = fe.getQ(fem, obs);
fem = fe.getQfull(fem);

fem = fe.FEMassemble(fem, 'output', 'matrices', 'verbose', true);

fem = fe.removeDirichlet(fem, 'sigmaBCL', sigmaBCL, 'sigmaBCR', sigmaBCR, ...
    'thicknessBCL', thkBCL, 'thicknessBCR', thkBCR, 'frequency', freq);

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
ylim([0.1 100]);
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
% for i = 1:fem.mesh.nt
%     r(:, i) = fem.mesh.Bk(:, :, i) * [1 / 3; 1 / 3] + fem.mesh.node(1:2, fem.elem2dofs(3, i));
% end
% DT = delaunay(r');
% X = tools.asColumn(r(1, :));
% Y = tools.asColumn(r(2, :));
% Z = tools.asColumn(abs(fem.Q.QJefull * sol.ue));
% L = logspace(ceil(log10(min(Z))), floor(log10(max(Z))), 50);

%%
figure(3);
plot.patchplotConst(fem.mesh.node, ...
    [fem.elem2dofs(1:3, :); ...
    tools.asRow(abs(fem.Q.QJefull * sol.ue))], ...
    @(x) x, 'none');
% hold on;
% plot.tricont(X, Y, DT, Z, L, 'w--');
% hold off;
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