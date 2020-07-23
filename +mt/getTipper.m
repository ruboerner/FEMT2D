function sol = getTipper(sol)
%getTipper Compute Tipper from Magnetic fields in E-polarization
%

if sol.epol
    assert(isfield(sol, 'Hz'));
    sol.T = sol.Hz ./ sol.Hy;
end