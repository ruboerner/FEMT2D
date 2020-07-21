function msh = getMesh(varargin)
%getMesh(varargin) read mesh file
% mesh = getMesh('filename', 'test.1', 'format', 'triangle');
%
p = inputParser;

addParameter(p, 'filename', '', @ischar);
addParameter(p, 'format', 'triangle', @ischar);
addParameter(p, 'scale', 1, @isnumeric);
addParameter(p, 'shift', 0, @isnumeric);
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

filename = p.Results.filename;
format = p.Results.format;
scale = p.Results.scale;
shift = p.Results.shift;
verbose = p.Results.verbose;

if verbose
    fprintf('\nMesh input:\n');
    fprintf('===========\n');
    fprintf('  Read mesh...');
    t0 = tic();
end
switch format
    case 'triangle'
        msh = mesh.getMeshFromTriangleFile(filename, scale, shift);        
    otherwise
        error('Not implemented.');
end
if verbose
    tRead= toc(t0);
    fprintf('done (%.2f s).\n', tRead);
end

msh = mesh.completeMesh(msh);

if verbose
    fprintf('\nMesh statistics:\n================\n');
    fprintf('  Triangles: %d\n  Edges: %d\n  Nodes: %d\n', ...
        msh.nt, msh.ne, msh.np);
    fprintf('  Boundary nodes: %d\n  Boundary edges: %d\n', ...
        length(msh.bdNodes), length(msh.bdEdges));
    fprintf('  Min triangle area: %.3e m^2\n  Max triangle area: %.3e m^2\n', ...
        min(msh.area), max(msh.area));
end
