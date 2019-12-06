function runall_setup(ifmissing)
% book_setup()
% book_setup('warn')
% book_setup('error')
% book_setup('ignore')
%
% Adds the necessary directories to Octave/Matlab's path and checks for
% Casadi and MPCTools. Argument specifies what to do if Casadi/MPCTools are not
% found.
narginchk(0, 1);
if nargin() < 1
    ifmissing = 'warn';
end

%% put paresto.m in the path
addpath('../..')

% Issue errors, warnings, etc.
switch ifmissing
case 'error'
    missing = @error;
case 'warn'
    missing = @warning;
case 'ignore'
    missing = @(varargin) [];
otherwise
    error('Invalid argument!');
end

%% Assume Casadi is in the path or print an error

if isempty(which('casadiMEX'))
    missing(['Casadi not found. Download from `casadi.org` and unzip ' ...
             'to a folder called ''casadi'' in your path.']);
end

end%function

