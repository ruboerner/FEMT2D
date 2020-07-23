function fem = removeDirichlet(fem, varargin)
%removeDirichlet remove Dirichlet values from system of linear equations
%
% fem = removeDirichlet(fem, varargin) removes Dirichlet condition from 
% system of linear equations
%
% Input:
% ======
%
% For a given frequency, removeDirichlet() computes electric and/or
% magnetic fields along the artificial boundaries of the computational
% domain. The conductivity structure is assumed to be a layered halfspace.
% However, the layers at both left and right boundaries may differ from
% each other.
%
% | 'sigmaBC', 'sigmaBCL', 'sigmaBCR' | conductivities in S/m at both, 
% |                                   | left or right boundaries
% | 'thicknessBC', thicknessBCL',     | layer thicknesses in m
% | 'thicknessBCR'                    |
% | 'frequency'                       | frequency in Hz
%
% (C) 2020 Ralph-Uwe Börner
%

p = inputParser;

addParameter(p, 'sigmaBC', [], @isnumeric);
addParameter(p, 'sigmaBCL', [], @isnumeric);
addParameter(p, 'sigmaBCR', [], @isnumeric);

addParameter(p, 'thicknessBC', [], @isnumeric);
addParameter(p, 'thicknessBCL', [], @isnumeric);
addParameter(p, 'thicknessBCR', [], @isnumeric);
addParameter(p, 'frequency', [], @isnumeric);

parse(p, varargin{:});

sigmaBC = p.Results.sigmaBC;
thicknessBC = p.Results.thicknessBC;

sigmaBCL = p.Results.sigmaBCL;
thicknessBCL = p.Results.thicknessBCL;
sigmaBCR = p.Results.sigmaBCR;
thicknessBCR = p.Results.thicknessBCR;
freq = p.Results.frequency;

assert(~isempty(sigmaBC) | ~isempty(sigmaBCL) | ~isempty(sigmaBCR), ...
    'Conductivity at boundary missing.');

if isempty(sigmaBC)
    
    if isempty(sigmaBCL)
        sigmaBCL = sigmaBCR;
        thicknessBCL = thicknessBCR;
    end
    
    if isempty(sigmaBCR)
        sigmaBCR = sigmaBCL;
        thicknessBCR = thicknessBCL;
    end
    
else
    
    [sigmaBCL, sigmaBCR, thicknessBCL, thicknessBCR] = deal(...
        sigmaBC, sigmaBC, thicknessBC, thicknessBC);
end



assert(length(sigmaBCL) == 1 + length(thicknessBCL), ...
    'Size of discrete boundary conductivities does not match number of thicknesses.');
assert(length(sigmaBCR) == 1 + length(thicknessBCR), ...
    'Size of discrete boundary conductivities does not match number of thicknesses.');


iw = 2i * pi * freq;
fem.app.frequency = freq;

fem.uDir = zeros(fem.ndofs, 1);

% For top and bottom boundaries, apply something like
%
% cs = spline([-1e4, 1e4], [0 -1 1 0]);
% xq = -1e4:100:1e4;
% plot(xq, ppval(cs, xq))

if strcmp(fem.application, 'MT')
    if strcmp(fem.polarization, 'epol') || strcmp(fem.polarization, 'both')
        fem.uDir = complex(zeros(fem.ndofs, 1));
        
        xleftbd = min(fem.coorddofs(1,fem.bddofs));
        xrightbd = max(fem.coorddofs(1,fem.bddofs));
        ztopbd = min(fem.coorddofs(2,fem.bddofs));
        zbotbd = max(fem.coorddofs(2,fem.bddofs));
        
        cs = spline([xleftbd xrightbd], ...
        [0 mt.getE1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, ztopbd) ...
        mt.getE1dMT(freq, 1 ./ sigmaBCR, thicknessBCR, ztopbd) ...
        0]);
    
        fem.uDir(fem.bddofstop) = ppval(cs, fem.coorddofs(1, fem.bddofstop));
        
        cs = spline([xleftbd xrightbd], ...
        [0 mt.getE1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, zbotbd) ...
        mt.getE1dMT(freq, 1 ./ sigmaBCR, thicknessBCR, zbotbd) ...
        0]);
    
        fem.uDir(fem.bddofsbot) = ppval(cs, fem.coorddofs(1, fem.bddofsbot));
        
        
        fem.uDir(fem.bddofsleft) = ...
            mt.getE1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, ...
            fem.coorddofs(2, fem.bddofsleft));
        fem.EnL = mt.getE1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, 0.0);
        
        fem.uDir(fem.bddofsright) = ...
            mt.getE1dMT(freq, 1 ./ sigmaBCR, thicknessBCR, ...
            fem.coorddofs(2, fem.bddofsright));
        
        fem.EnR = mt.getE1dMT(freq, 1 ./ sigmaBCR, thicknessBCR, 0.0);
        
        fem.uDir(fem.bddofstop) = ...
            mt.getE1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, ...
            fem.coorddofs(2, fem.bddofstop));
        fem.uDir(fem.bddofsbot) = ...
            mt.getE1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, ...
            fem.coorddofs(2, fem.bddofsbot));
        
        fem.Ac = fem.Null' * (fem.Ke + iw * fem.Me) * fem.Null;
        fem.fc = fem.Null' * (fem.uDir - (fem.Ke + iw * fem.Me) * fem.uDir);
    end
    
    if strcmp(fem.polarization, 'hpol') || strcmp(fem.polarization, 'both')
        fem.uDirh = complex(zeros(fem.ndofs, 1));
        fem.uDirh(fem.stripdofs) = ...
            mt.getH1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, ...
            fem.coorddofs(2, fem.stripdofs));
        
        fem.uDirh(fem.bddofsleft) = ...
            mt.getH1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, ...
            fem.coorddofs(2, fem.bddofsleft));
        fem.uDirh(fem.bddofsright) = ...
            mt.getH1dMT(freq, 1 ./ sigmaBCR, thicknessBCR, ...
            fem.coorddofs(2, fem.bddofsright));
        
        cs = spline([xleftbd xrightbd], ...
        [0 mt.getH1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, zbotbd) ...
        mt.getH1dMT(freq, 1 ./ sigmaBCR, thicknessBCR, zbotbd) ...
        0]);
    
        fem.uDirh(fem.bddofsbot) = ppval(cs, fem.coorddofs(1, fem.bddofsbot));
        
        fem.HnL = mt.getH1dMT(freq, 1 ./ sigmaBCL, thicknessBCL, 0.0);
        fem.HnR = mt.getH1dMT(freq, 1 ./ sigmaBCR, thicknessBCR, 0.0);
        
        fem.Ach = fem.NullStrip' * (fem.Kh + iw * fem.Mh) * fem.NullStrip;
        fem.fch = fem.NullStrip' * (fem.uDirh - (fem.Kh + iw * fem.Mh) * fem.uDirh);
    end
end
