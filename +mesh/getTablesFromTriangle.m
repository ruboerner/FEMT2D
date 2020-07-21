function [p, t, v, marker] = getTablesFromTriangle(filename)
%getTablesFromTriangle(filename)
%
nodeFile = [filename '.node'];
elemFile = [filename '.ele' ];

fid = fopen(nodeFile, 'r');

slurp = fscanf(fid, '%i', [1 4]);
nNodes = slurp(1);
ncol = 3 + slurp(3) + slurp(4);

p = zeros(ncol - 1, nNodes);
slurp = fscanf(fid, '%f', [ncol nNodes]);
p = slurp(2:ncol, :);

fclose(fid);

fid = fopen(elemFile, 'r');

slurp = fscanf(fid, '%d', [1 3]);
nElem = slurp(1);
ncol = 4 + slurp(3);

% T = zeros(4, nElem);

slurp = fscanf(fid, '%f', [ncol nElem]);
T_ = slurp(1:ncol, :);
    
t = T_(2:4, :);
if T_(1,1) == 0
    t = t + 1;
end
v = T_(5, :);
fclose(fid);

marker = p(3, :);
p = p(1:2, :);

