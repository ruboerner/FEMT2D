function sol = postProcessing(fem, sol)
%postProcessing Helper function to capsule all post-processing steps into one function
%

sol = mt.getFields(fem, sol);
sol = mt.getTipper(sol);
sol = mt.getImpedance(sol);
sol = mt.getRhoa(fem, sol);
sol = mt.getPhase(sol);
