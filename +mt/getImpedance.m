function sol = getImpedance(sol)
%getImpedance Compute impedance from electric and magnetic fields
%

if sol.epol
    assert(isfield(sol, 'Ex') && isfield(sol, 'Hy'));
    sol.Zxy = sol.Ex ./ sol.Hy;
end

if sol.hpol
    assert(isfield(sol, 'Ey') && isfield(sol, 'Hx'));
    sol.Zyx = sol.Ey ./ sol.Hx;
end