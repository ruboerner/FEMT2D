function patchplotU(fem, u, col)
% Plotfunktion zur Darstellung der Loesung. Jedes Dreieck der
% Triangulierung wird gezeichnet als 3D 'Patch' durch die MATLAB-Funktion
% patch. Koordinaten der Ecken werden aus dem p-Feld ausgelesen. Zum
% Zeichnen werden alle Dreieck abgelaufen (t-Feld).
%
% Oliver Ernst, 29.05.2009
% Numerische Simulation mit finiten Elementen
% Sommersemester 2009
% TU Bergakademie Freiberg
% Institut fuer Numerische Mathematik und Optimierung

% clf
if nargin < 4
    col = 'none';
end

trisurf(fem.mesh.tri2node(1:3,:)', fem.mesh.node(1,:)', fem.mesh.node(2,:)', u, 'EdgeColor', col)

view(2);
axis equal tight ij
box on
