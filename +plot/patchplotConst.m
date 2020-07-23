function patchplotConst(p, t, fun, color)
%patchPlotConst Plot coefficients on a triangulation
%
% patchplotConst(p, t, fun, color)
%
% Plots a triangualar mesh defined by np points giben in the array p(1:2, 1:np).
% The triangular mesh is further defined by the T(1:3, ne) coefficients
% The 4-th row of T, i.e, T(4, 1:ne), contains the quantity to be visualized.
% T(4, 1:ne) can be transformed by the function fun.
% The triangular grid can be colorized by defining a color string.
%


if nargin < 3
    fun = @log10;
end

if nargin < 4
    color = 'none';
end

% clf
patch( ...
    reshape(p(1, t(1:3, :)), 3, []), ...
    reshape(p(2, t(1:3, :)), 3, []), ...
    fun(t(end, :)), 'EdgeColor', color, 'Linewidth', 0.01);
% axis equal tight ij
% box on
