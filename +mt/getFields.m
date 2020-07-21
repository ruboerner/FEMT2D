function sol = getFields(fem, sol)

both = strcmp(fem.polarization, 'both');

if strcmp(fem.polarization, 'epol') || both
    iw = -2i * pi * fem.app.frequency;
    Ex = fem.Q.Qe   * sol.ue;
    Hy = fem.Q.QeHy * sol.ue / iw;
    Hz = fem.Q.QeHz * sol.ue / iw;
    sol.epol = true;
    sol.Ex = Ex;
    sol.Hy = Hy;
    sol.Hz = Hz;
end

if strcmp(fem.polarization, 'hpol') || both
    Ey = fem.Q.QhEy * sol.uh;
    Ez = fem.Q.QhEz * sol.uh;
    Hx = fem.Q.Qh   * sol.uh;
    sol.hpol = true;
    sol.Hx = Hx;
    sol.Ey = Ey;
    sol.Ez = Ez;
end