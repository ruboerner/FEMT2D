function plotMT(varargin)

p = inputParser();

addParameter(p, 'fem', struct([]), @isstruct);
addParameter(p, 'mesh', struct([]), @isstruct);
addParameter(p, 'sol', struct([]), @isstruct);

addParameter(p, 'obs', [], @isnumeric);
addParameter(p, 'periods', [], @isnumeric);

addParameter(p, 'profile', [], @ischar);
addParameter(p, 'section', [], @ischar);
addParameter(p, 'sounding', [], @ischar);
addParameter(p, 'station', [], @isnumeric);

% profile: rhoa+phase, rhoa, phase, tipper
%
% section: conductivites, current density abs(Jx), fields
%
% sounding: rhoa+phase vs. period
%

parse(p, varargin{:});

f = @(x) getfield(p.Results, x); %#ok

fem = f('fem');
mesh = f('mesh');
sol = f('sol');
obs = f('obs');
periods = f('periods');


profile = f('profile');
section = f('section');
sounding = f('sounding');
site = f('station');

if strcmp(profile, 'rhoa+phase')
    title('Profile of rhoa and phase');
    xobs = obs(1, :);
    subplot(2, 1, 1);
    semilogy(xobs, sol.rhoaxy, '.-', xobs, sol.rhoayx, '.-');
    
    mi = floor(log10(min(min(sol.rhoaxy, sol.rhoayx))));
    ma = ceil(log10(max(max(sol.rhoaxy, sol.rhoayx))));
    
    ylim(10.0 .^ [mi ma]);
    grid();
    xlabel('y in km');
    title(['Apparent resistivities ' sprintf(...
        'T = %.1e s', 1.0 ./ fem.app.frequency)]);
    legend('E-pol', 'H-Pol');
    subplot(2, 1, 2);
    plot(xobs, sol.phixy, '.-', xobs, sol.phiyx, '.-');
    ylim([0 90]);
    grid();
    xlabel('y in km');
    title(['Impedance phases ' sprintf(...
        'T = %.1e s', 1.0 ./ fem.app.frequency)]);
    legend('E-pol', 'H-Pol');
end

if strcmp(section, 'Jx')
    plot.patchplotConst(fem.mesh.node, ...
        [fem.elem2dofs(1:3, :); ...
        tools.asRow(abs(fem.Q.QJefull * sol.ue))], ...
        @(x) x, 'none');
    title('Current density |Jx|');
    colormap(gca, jet(256));
    h = colorbar();
    ylabel(h, '|J| in A/m^2');
    axis equal tight ij
end

if strcmp(profile, 'tipper')
    xobs = obs(1, :);
    plot(xobs, real(sol.T), 'r--', xobs, imag(sol.T), 'b--');
    grid();
    title(['Tipper ' sprintf('T = %.1e s', ...
        1.0 ./ fem.app.frequency)]);
    xlabel('y in km');
    ylabel('magn. transfer function');
    legend('Re', 'Im');
end

if strcmp(sounding, 'rhoa+phase')
    assert(length(site) == 1, 'plotMT: Sounding station required');
    subplot(1, 2, 1);
    rhoaxy = sol.rhoaxy(site, :);
    rhoayx = sol.rhoayx(site, :);
    
    loglog(periods, ...
        rhoaxy, 'r.-', ...
        periods, ...
        rhoayx, 'b.-');
    ylim(10.^[floor(log10(min(min(rhoaxy, rhoayx)))), ...
        ceil(log10(max(max(rhoaxy, rhoayx)))) ]);
    xlabel('T in s');
    ylabel('rhoa in Ohm*m');
    grid();
    legend('E-Pol.', 'H-Pol.');
    subplot(1, 2, 2);
    phixy = sol.phixy(site, :);
    phiyx = sol.phiyx(site, :);
    
    semilogx(periods, ...
        phixy, 'r.-', ...
        periods, ...
        phiyx, 'b.-');
    xlabel('T in s');
    ylabel('phase in degrees');
    grid();
    ylim([0 90]);
    legend('E-Pol.', 'H-Pol.');
end

if strncmp(section, 'cond', 4)
    plot.patchplotConst(fem.mesh.node, ...
        [fem.mesh.tri2node; tools.asRow(fem.sigma)], @log10, [0.9 0.9 0.9]);
    xlabel('y in m');
    ylabel('z in m');
    title('Conductivity in S/m');
    colormap(gca, jet(256));
    h = colorbar();
    axis equal tight ij
    ylabel(h, 'log_{10}(\sigma/({S/m}))');
end

if strncmp(section, 'subd', 4)
    plot.plotSubdomains(mesh, 'meshcolor', 'none');
    title('Subdomain order');
    xlabel('y in m');
    ylabel('z in m');
    axis equal tight ij
end