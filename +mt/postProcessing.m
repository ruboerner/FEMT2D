function sol = postProcessing(fem, sol)

sol = mt.getFields(fem, sol);
sol = mt.getTipper(sol);
sol = mt.getImpedance(sol);
sol = mt.getRhoa(fem, sol);
sol = mt.getPhase(sol);
