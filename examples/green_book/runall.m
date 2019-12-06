function failures = runall(doplots)
% failures = runall([showplots='off'])
%
% Runs all of the figure scripts and returns a cell array of errors for any
% scripts that failed.
%
% If showplots is true, plots will be displayed as they are created (note that
% they will steal focus, which can be annoying). Otherwise, they will remain
% hidden until the script finishes (at which point they will also steal focus).
% Finally, if showplots is the string 'off', then they will never be shown.
narginchk(0, 1);
clearplots = false();
if nargin() < 1
    doplots = 'off';
end
if isequal(doplots, 'off')
    doplots = false();
    clearplots = true();
end

% Make sure user has the required packages.
runall_setup('error');

% Set some local printing options if octave.
if isOctave()
    more('off');
    page_output_immediately(true(), 'local');
    page_screen_output(false(), 'local');
end

% Start with a dummy NLP to display the IPOPT splash message.
x = casadi.SX.sym('x');
nlp = struct('x', x, 'f', x.^2);
ipopt = casadi.nlpsol('ipopt', 'ipopt', nlp, ...
                      struct('print_time', false(), ...
                             'ipopt', struct('print_level', 0)));
ipopt('x0', 0);

% Switch off plotting.
DefaultFigureVisible = get(0, 'DefaultFigureVisible');
if ~doplots
    set(0, 'DefaultFigureVisible', 'off'); % Don't show figures.
end
initialfigures = get(0, 'children');
getnewplots = @() setdiff(get(0, 'children'), initialfigures);

% Now actually run the examples.
fprintf('\nRunning book scripts.\n\n');
bookfigures = {
	       %% Only runs unique figure scripts.
	       'adsone';
	       'adsall';
	       'batch_data_solution';
	       'bvsm';
	       'bvsm_red';
	       'estdiff';
	       'fitrtd';
	       'hbv_det';
	       'hbv_red';
	       'react2rev';
	       'ABC';
	       'Sfirstorder';
};
failures = cell(0, 1);
for i = 1:length(bookfigures)
    fig = bookfigures{i};
    fprintf('* %s ... ', fig);
    [err, figtime] = runmain(fig);
    if isempty(err)
        fprintf('success');
    else
        fprintf('FAILURE');
        failures{length(failures) + 1} = struct('figure', fig, 'err', err);
    end
    fprintf(' (%.4g s)\n', figtime);
    if clearplots
        close(getnewplots()); % Close figures to save memory.
    end
end

% Make figures show up.
if ~doplots && ~clearplots
    newplots = getnewplots();
    for i = 1:length(newplots)
        set(newplots(i), 'Visible', 'on');
    end
end

% Switch plots back on.
set(0, 'DefaultFigureVisible', DefaultFigureVisible);

% Display failures if output not requested.
if nargout() == 0
    for i = 1:length(failures)
        f = failures{i};
        fprintf('*** %s ***\n', f.figure);
        disp(f.err);
        
        % Print the stack trace, but skip the last two (they will always be
        % runall, runmain, and runscript_octave/matlab.
        for j = 1:(length(f.err.stack) - 3)
            disp(f.err.stack(j));
        end
    end
end

end%function

function [err, figtime] = runmain(fig)
    % [err, figtime] = runmain(fig)
    %
    % Changes to the directory of of the given figure and runs main.m, returning
    % any error message. Afterwards, changes back to the original directory.
    ## olddir = cd(fig);
    ## dircleanup = onCleanup(@() cd(olddir));
    ## oldpath = path();
    ## pathcleanup = onCleanup(@() path(oldpath));
    
    figure(); % Avoids issues with reusing old figures.
    figtime = tic();
    if isOctave()
      pkg load statistics
      err = runscript_octave(fig);
    else
      err = runscript_matlab(fig);
    end
    figtime = toc(figtime);
end%function

function err = runscript_matlab(script)
    % err = runscript_matlab(script)
    %
    % Attempts to run a script without displaying anything and returns any error
    % message.
    try
        evalc(script);
        err = [];
    catch err
        % Pass.
    end
end%function

function err = runscript_octave(script)
    % err = runscript_octave(script)
    %
    % Attempts to run a script and returns any error message.
    %
    % Most output should be suppressed, although we can't hit everything.
    clear('-x', 'script'); % Clear any lingering function definitions.
    try
        evalc(script);
        err = [];
    catch err
        % Pass.
    end
end%function

function noop(varargin)
    % Does nothing.
end%function

