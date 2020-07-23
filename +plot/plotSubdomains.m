function plotSubdomains(mesh, varargin)
%plotSubdomains(mesh)
%
% (C) 2020 Ralph-Uwe Börner

p = inputParser;

addParameter(p, 'meshcolor', 'none', @ischar);
parse(p, varargin{:});

meshcolor = p.Results.meshcolor;

subdomains = unique(mesh.tri2subdomain);
nsubdomains = length(subdomains);
smallestsubdomain = min(subdomains);

plot.patchplotConst(mesh.node, [mesh.tri2node; mesh.tri2subdomain], @(x) x, meshcolor);
set(gca, 'CLim', -0.5 + [smallestsubdomain, nsubdomains+smallestsubdomain]);
colormap(gca, jet(nsubdomains));
h = colorbar();
ylabel(h, 'subdomain number');
set(h, 'YTick', unique(mesh.tri2subdomain));
