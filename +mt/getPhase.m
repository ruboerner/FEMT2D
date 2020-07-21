function sol = getPhase(sol)

if sol.epol
    assert(isfield(sol, 'Zxy'));
    sol.phixy = atan(imag(sol.Zxy) ./ real(sol.Zxy)) * 180 / pi;
end
if sol.hpol
    assert(isfield(sol, 'Zyx'));
    sol.phiyx = atan(imag(sol.Zyx) ./ real(sol.Zyx)) * 180 / pi;
end
