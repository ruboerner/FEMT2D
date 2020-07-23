function idxDOFS = finddofsinsubdomain(fem, isubdomain)
%finddofsinsubdomain() find DOFS associated with a subdomain
%
% indexDOFS = finddofsinsubdomain(fem, subdomainnumber)
%
% Given a mesh stored in fem.mesh, finddofsinsubdomain() collects all
% indices of dofs which correspond to a subdomain indicated by
% subdomainnumber
%
% Input:
% ======
% fem - fem structure, see FEMproblem()
% subdomainnumber - number of subdomain for which the DOFS are desired
%
% (C) 2020 Ralph-Uwe Börner
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

