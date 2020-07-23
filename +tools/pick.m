function element = pick(which, varargin)
%pick Picks one of the given arguments.
%
% SYNTAX
%
%   X = pick(index, arg1[, arg2[, ...]])
%
% INPUT/OUTPUT PARAMETERS
%
%   X     ... Selected argument.
%   index ... Index of the argument that should be returned. E.g.
%             specify '2' to assign arg2 to X.
%   arg1  ... See description of 'index'.
%   arg2  ... See description of 'index'.
%
% (C) Martin Afanasjew
%

assert(which >= 1 && which <= length(varargin), ...
    'Invalid index specified. You can only pick from the arguments you have supplied.');

element = varargin{which};
end
