clearvars();
close('all');
%clc();

%%
meshFile = 'meshes/commemi2d4.1';

mesh = mesh.getMesh('filename', meshFile, 'format', 'triangle', ...
    'shift', 0, 'scale', 1000, 'verbose', true);

%%
plot.plotSubdomains(mesh, 'meshcolor', 'none');
xlabel('y in m');
ylabel('z in m');
axis equal tight ij

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

En = fem.EnL;
Hn = fem.HnL;

%%
figure(2);
subplot(2, 1, 1);
semilogy(xobs, rhoaxy, '.-', xobs, rhoayx, '.-');
ylim([0.1 1000]);
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

%% COMMEMI comparison
% Average(1) values are used here for
% T = 100 s.
%
x = [-10, -7, -6, -5, 2, 3.5, 5] * 1000;
% rhoaE = [12.69, 11.7, 8.68, 6.76, 6.49, 6.11, 6.18];
rhoaE = [37.2, 28.28, 23.44, 19.19, 18.02, 20.81, 27.7];
% rhoaH = [11.41, 11.30, 8.97, 6.78, 6.85, 6.33, 5.84];
rhoaH = [74.24, 66.69, 32.58, 7.11, 4.39, 8.52, 39.06];

plot(xobs, rhoaxy, 'r.-', x, rhoaE, 'ro-', ...
    xobs, rhoayx, 'b.-', x, rhoaH, 'bo-');
xlim([-12., 7] * 1000);
grid();

%%
% E-Pol. 
xe = [-10, -7, -6, -5, 0, 2, 3.5, 5, 8, 16] * 1000;

Exr = [0.846, 0.823, 0.82, 0.808, 0.801, 0.806, 0.8, 0.824, 0.841, 0.869];
Exi = [0.181,0.220,0.233,0.252,0.279,0.267,0.252,0.230,0.197,0.149];
plot(xobs, imag(Ex ./ fem.EnL), ...
    xe, -Exi, 'o')
xlim(1000 * [-12.0, 18.0]);
%%
Hyr = [1.049, 1.202, 1.315, 1.428, 1.534, 1.462, 1.370, 1.208, 1.076, 1.006];
Hyi = [0.010,-0.079,-0.179,-0.245,-0.314,-0.273,-0.201,-0.102,-0.012,0.018];

plot(xobs, imag(Hy ./ fem.HnL), ...
    xe, -Hyi, 'o')
xlim(1000 * [-12.0, 18.0]);

%%
Hzr = [-0.287, -0.396, -0.399, -0.350, 0.077, 0.219, 0.318, 0.326, 0.240, 0.116];
Hzi = [0.184,0.261,0.275,0.227,-0.045,-0.135,-0.198,-0.214,-0.149,-0.054];
plot(xobs, imag(Hz ./ fem.HnL), ...
    xe, -Hzi, 'o')
xlim(1000 * [-12.0, 18.0]);

%% 
% H-Pol.
Eyr = [1.131,1.066,0.734,0.343,0.227,0.277,0.394,0.789,0.893,0.908];
plot(xobs, real(Ey ./ fem.EnL), ...
    xe, -Eyr, 'o')
xlim(1000 * [-12.0, 18.0]);

%%
Eyi = [0.054,0.046,0.026,-0.004,-0.017,-0.008,0.006,0.042,0.053,0.054];
plot(xobs, imag(Ey ./ fem.EnL), ...
    xe, Eyi, 'o')
xlim(1000 * [-12.0, 18.0]);


