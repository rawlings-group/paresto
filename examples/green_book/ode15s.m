function [t_out, x] = ode15s (f, t, x0, options)

  if (nargin == 3 || nargin == 4)

    global __ode15s_rhs_fcn__;
    global __ode15s_event_fcn__;

    __ode15s_rhs_fcn__ = f;

    use_dasrt = false;
    if (nargin == 4)
      if (isfield (options, "Events") && ! isempty (options.Events))
	use_dasrt = true;
	__ode15s_event_fcn__ = options.Events;
	options_fcn = @dasrt_options;
      else
	options_fcn = @lsode_options;
      endif
      if (isfield (options, "RelTol") && ! isempty (options.RelTol))
	options_fcn ("relative tolerance", options.RelTol)
      endif
      if (isfield (options, "AbsTol") && ! isempty (options.AbsTol))
	options_fcn ("absolute tolerance", options.AbsTol)
      endif
    endif

    if (use_dasrt)
      xdot0 = f (t(1), x0);
      [x, xdot, t_out] = dasrt ("__dasrt_rhs__", "__dasrt_event_fcn__",
				x0, xdot0, t);
    else
      x = lsode ("__lsode_rhs__", x0, t);
      t_out = t;
    endif

  else
    usage ("[t_out, x] = ode15s (f, t, x0, options)");
  endif

endfunction

% ----------------------------------------------------------------------

function xdot = __lsode_rhs__ (x, t)

  global __ode15s_rhs_fcn__;

  xdot = __ode15s_rhs_fcn__ (t, x);

endfunction

% ----------------------------------------------------------------------

function resid = __dasrt_rhs__ (x, xdot, t)

  global __ode15s_rhs_fcn__;

  resid = xdot - __ode15s_rhs_fcn__ (t, x);

endfunction

% ----------------------------------------------------------------------

function retval = __dasrt_event_fcn__ (x, t)

  global __ode15s_event_fcn__;

  retval = __ode15s_event_fcn__ (t, x);

endfunction
