function idxDOFS = finddofsinsubdomain(fem, isubdomain)
%finddofsinsubdomain(fem, isubdomain) find dofs indices associated with a
%subdomain
%

ns = unique(fem.mesh.tri2subdomain);
assert((isubdomain >= min(ns)) & (isubdomain <= max(ns)), ...
    'Subdomain number out of range.');

isub = find(fem.mesh.tri2subdomain == isubdomain);
idxNodes = unique(reshape(fem.mesh.tri2node(1:3, isub), [], 1));
idxEdges = fem.mesh.np + unique(reshape(fem.mesh.tri2edge(1:3, isub), [], 1));

idxDOFS = idxNodes;
if fem.order == 2
    idxDOFS = [idxDOFS; idxEdges];
end

