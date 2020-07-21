function patchplotConst(p, t, fun, color)
%function patchplotConst(p, t, fun, color)
%
% Plotfunktion zur Darstellung stueckweise konstanter Koeffizienten.
% Jedes Dreieck der Triangulierung wird gezeichnet als 3D 'Patch' durch
% die MATLAB-Funktion patch.
% Koordinaten der Ecken werden aus dem p-Feld ausgelesen.
% Der Koeffizient befindet sich in der letzten Zeile der Elementtabelle t.
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
