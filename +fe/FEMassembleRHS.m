function f = FEMassembleRHS(fem, R)
%FEMassembleRHS assemble RHS of linear system for given source
%

% Identify index and barycentric coordinates of that element
% to which source coordinate belongs
%
[indEl, bary] = tsearchn(fem.mesh.node(1:2, :)', ...
    fem.mesh.tri2node(1:3, :)', ...
    [asColumn(R(1, :)), asColumn(R(2, :))]);
idx = fem.elem2dofs(1:3, indEl);
f = zeros(fem.ndofs, 1);
f(idx) = R(1) * bary;
