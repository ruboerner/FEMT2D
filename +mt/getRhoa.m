function sol = getRhoa(fem, sol)

omega = 2 * pi * fem.app.frequency;
mu0 = pi * 4e-7;

if sol.epol
    assert(isfield(sol, 'Zxy'));
    sol.rhoaxy = abs(sol.Zxy).^2 / mu0 / omega;
end
if sol.hpol
    assert(isfield(sol, 'Zyx'));
    sol.rhoayx = abs(sol.Zyx).^2 / mu0 / omega;
end
