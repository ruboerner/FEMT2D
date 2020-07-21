function sol = getTipper(sol)

if sol.epol
    assert(isfield(sol, 'Hz'));
    sol.T = sol.Hz ./ sol.Hy;
end