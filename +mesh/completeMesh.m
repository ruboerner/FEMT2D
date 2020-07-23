function mesh = completeMesh(mesh)
%completeMesh calculates missing quantities from a triangular mesh
%
% completeMesh computes point and edge tables as well as affine maps
% from the information stored in the mesh struct as provided by, e.g.,
% mesh generators like Triangle.
%
% Input:
% =====
% A triangular mesh consisting of a triangle-to-node map stored in the
% struct mesh
%
% Output:
% =======
%
% completeMesh augments the struct mesh.
%
% (C) in its original form by Martin Afanasjew
%

% The following array associates the three edges of a triangle with the
% three local vertex indices (see description above).
%
edgesByVert = [ 2 3;  % edge 1
    3 1;  % edge 2
    1 2]; % edge 3

% Replicate the rows of the element table according to the edgesByVert
% table. This yields an array the columns of which contain the global
% vertex indices making up each of the three edges (6 in all).
%
edges = mesh.tri2node(edgesByVert', :);

% Reshaping this array as 2 x 3*nTri results in an array whose columns
% contain the 2 global vertex indices (in order) of all edges in the mesh,
% where interior edges occure twice.
%
edges = reshape(edges, 2, 3 * size(mesh.tri2node, 2));

% Sort this array so that all columns are ascending giving unique ordering
% to each edge.
%
edges = sort(edges, 1);

% Now we use Matlab's 'unique' command to remove multiple occurrences of
% edges and obtain the updated index to each single occurrence. 'Unique' can
% work on rows, so we transpose first.
%
[edges, ~, idx] = unique(edges', 'rows');

% Now transpose 'edges' to get one edge in each column and group the
% entries of 'idx' three for each element.
%
mesh.edge2node = edges';
mesh.tri2edge = reshape(idx, 3, size(mesh.tri2node, 2));

% Calculate number of triangles, nodes and edges
%
mesh.nt = size(mesh.tri2node, 2);
mesh.np = size(mesh.node, 2);
mesh.ne = size(mesh.edge2node, 2);

elem = mesh.tri2node';

totalEdge = [elem(:,[2, 3]); elem(:,[3, 1]); elem(:,[1, 2])];

totalEdge = sort(totalEdge, 2);
[i, j, s] = find(...
    sparse(double(totalEdge(:,1)), double(totalEdge(:,2)), 1));
bdEdge = [i(s == 1 ), j(s == 1)];
isBdNode = false(max(mesh.tri2node(:)), 1);
isBdNode(bdEdge(:)) = true;
bdNode = find(isBdNode);

[~, ~, n] = unique(totalEdge, 'rows');
counts = accumarray(n(:), 1);

mesh.bdEdges = find(counts == 1); %
% mesh.bdEdges = find(all(ismember(mesh.edge2node, bdEdge'), 1))';
mesh.bdNodes = bdNode;

get_Bk = @(i) [ ...
    mesh.node(1, mesh.tri2node(1, i)) - mesh.node(1, mesh.tri2node(3, i)) mesh.node(1, mesh.tri2node(2, i)) - mesh.node(1, mesh.tri2node(3, i)); ...
    mesh.node(2, mesh.tri2node(1, i)) - mesh.node(2, mesh.tri2node(3, i)) mesh.node(2, mesh.tri2node(2, i)) - mesh.node(2, mesh.tri2node(3, i)) ...
    ];

Bk = zeros(2, 2, mesh.nt);
Binv = zeros(2, 2, mesh.nt);
detBk = zeros(mesh.nt, 1);
area = zeros(mesh.nt, 1);

for ii = 1:mesh.nt
    Bk(:, :, ii) = get_Bk(ii);
    Binv(:, :, ii) = inv(Bk(:, :, ii)');
    detBk(ii) = det(Bk(:, :, ii));
    if detBk(ii) < 0
        mesh.tri2node([2, 3], ii) = mesh.tri2node([3, 2], ii);
    end
    area(ii) = 0.5 * abs(detBk(ii));
end

mesh.Bk = Bk;
mesh.Binv = Binv;
mesh.detBk = detBk;
mesh.area = area;
