classdef paresto < handle
  properties
    % Print level (0 or missing field; no output)
    print_level
    % Method: Direct collocation or multiple shooting
    transcription
    % NLP solver plugin
    nlp_solver
    % NLP solver options
    nlp_solver_options
    % NLP presolver plugin
    nlp_presolver
    % NLP presolver options
    nlp_presolver_options
    % NLP solver integrator options (only used for the sundials transcription
    % options that uses sundials as the integrator during optimization).
    nlp_solver_integrator_options
    % Dynamic model
    model
    % ODE/DAE simulator
    simulator
    % Number of experiments
    nsets
    % Total number of data points in all experiments
    ndata
    % Number of measurement times - 1
    N
    % Number of states
    nx
    % Number of parameters
    np
    % Number of algebraic variables
    nz
    % Number of measurements or setpoints
    nd
    % Discrete-time dynamics
    daefun
    % Stage function
    stagefun
    % Output function
    outfun
    % Parameter function
    parfun
    % Collocation equations
    dynfun
    % Degree of interpolating polynomial
    ord
    % Roots of the collocation polynomial
    tau_root
    % Number of NLP decision variables
    nw
    % Components in w being estimated
    thetaind
    % All time points with measurements
    tout
    % Mapping to NLP decision vector
    to_w
    % Mapping from NLP decision vector
    from_w
    % NLP solver instance
    solver
    % NLP presolver instance
    presolver
    % Calculate parametric sensitivities of solver
    fsolver
  end
  methods
    function self = paresto(model)

      % Constructor

      % Control output to screen; default to no diagnostic output
      % Log message with timings
      if ~isfield(model, 'print_level')
        model.print_level = 1;
      end
      if (model.print_level > 0)
        msg = @(m) fprintf('paresto.paresto (t=%g ms): %s\n', 1000*toc, m);
        tic;
      else
	msg = @(m) fprintf('');
	model.nlp_solver_options.ipopt.print_level = 0;
	model.nlp_solver_options.print_time = false;
      end
      self.print_level = model.print_level;
      % Fields are empty by default
      f = {'x', 'z', 'p', 'y'};
      for i=1:numel(f)
        if ~isfield(model, f{i})
          model.(f{i}) = {};
        end
      end
      self.model = model;

      % Number of data sets (1 by default)
      if isfield(model, 'nsets')
        self.nsets = model.nsets;
      else
        self.nsets = 1;
      end

      % Total number of measurement points in all data sets
      if isfield(model, 'ndata')
        self.ndata = model.ndata;
      end

      % Order of interpolating polynomials
      if isfield(model, 'ord')
        self.ord = model.ord;
      else
        self.ord = 3;
      end

      % Get method
      msg('NLP transcription');
      if isfield(model, 'transcription')
        self.transcription = model.transcription;
      else
        self.transcription = 'simultaneous';
      end

      % NLP presolver
      if isfield(model, 'nlp_presolver')
        self.nlp_presolver = model.nlp_presolver;
      else
        self.nlp_presolver = [];
      end

      % NLP presolver options
      if isfield(model, 'nlp_presolver_options')
        self.nlp_presolver_options = model.nlp_presolver_options;
      else
        self.nlp_presolver_options = struct;
      end

      % NLP solver
      if isfield(model, 'nlp_solver')
        self.nlp_solver = model.nlp_solver;
      else
        self.nlp_solver = 'ipopt';
      end

      % NLP solver options
      if isfield(model, 'nlp_solver_options')
        self.nlp_solver_options = model.nlp_solver_options;
      else
        self.nlp_solver_options = struct;
      end

      % NLP solver integrator options.
      if isfield(model, 'nlp_solver_integrator_options')
        self.nlp_solver_integrator_options = model.nlp_solver_integrator_options;
      else
        self.nlp_solver_integrator_options = struct;
      end

      % Have NLP base class calculate multipliers
      self.nlp_solver_options.calc_lam_x = true;
      self.nlp_solver_options.calc_lam_p = true;

      % Do bound projection
      self.nlp_solver_options.bound_consistency = true;

      % Construct symbolic expressions for the dynamic model
      msg('DAE modeling');
      self.modeling();

      % Collocation equations
      msg('Collocation equations');
      self.collocation()

      % NLP transcription
      self.transcribe()

      % Prepare the sensitivity analysis
      msg('NLP sensitivity equations');
      self.fsolver = self.solver.forward(numel(self.thetaind));
      % Done intitializing
      msg('Initialization complete');
    end

    function retval = modeling(self)
      % Modeling stage: Create CasADi expressions from user data
      model = self.model;

      % Time
      t = casadi.SX.sym('t');

      % Free parameters
      pp = struct;
      [p,pp] = self.str2sym('p', pp);

      % Continuous variables
      yy = struct;
      [x,yy] = self.str2sym('x', yy); % Differential states
      [z,yy] = self.str2sym('z', yy); % Algebraic variables
      [d,yy] = self.str2sym('d', yy); % Measurements, data

      % Additional outputs
      y_def = self.fun2sym('h', t, yy, pp);
      [~,yy] = self.str2sym('y', yy);
      assert(numel(model.y)==numel(y_def));
      for i=1:numel(model.y)
        yy.(model.y{i}) = y_def(i);
      end

      % DAE
      ode = self.fun2sym('ode', t, yy, pp);
      alg = self.fun2sym('alg', t, yy, pp);
      % Check dg/dz
      if size(alg, 1)>0
        dalgdz = jacobian(alg, z);
        determinant_dalgdz = det(dalgdz);
        if is_zero(determinant_dalgdz)
          warning('The Jacobian of the algebraic equations wrt z is singular. You may have formulated a high index DAE.');
        end
      end

      % Least squares objective
      lsq = self.fun2sym('lsq', t, yy, pp);

      % Dimensions
      self.nx = numel(x);
      self.nz = numel(z);
      self.np = numel(p);
      self.nd = numel(d);

      % Create a simulator
      self.tout = model.tout(:)'; % Force row vector.
      self.N = numel(self.tout)-1;
      if self.nx > 0
        dae = struct('t', t, 'x', x, 'p', [p; d], 'z', z, 'ode', ode, 'alg', alg);
	grid = self.tout;
	t0 = grid(1);
        % Use CVODES for ODEs, IDAS for DAEs
        if self.nz==0
          plugin = 'cvodes';
        else
          plugin = 'idas';
        end
        self.simulator = casadi.integrator('simulator', plugin, dae, t0, grid);
      end

      % Function evaluated at each collocation point
      self.daefun = casadi.Function('daefun', {t, x, z, p, d}, {ode, alg}, ...
                                    {'t', 'x', 'z', 'p', 'd'}, {'ode', 'alg'});

      % Function evaluated at each measurement
      self.stagefun = casadi.Function('stagefun', {t, x, z, p, d}, {lsq, alg}, ...
                              {'t', 'x', 'z', 'p', 'd'}, {'lsq', 'alg'});

      % Output function
      self.outfun = casadi.Function('outfun', {t, x, z, p, d}, struct2cell(yy),...
                                    {'t', 'x', 'z', 'p', 'd'}, fieldnames(yy));

      % Parameter function
      self.parfun = casadi.Function('parfun', {p}, struct2cell(pp),...
                                    {'p'}, fieldnames(pp));
    end

    function [x, z] = simulate(self, d, xin, pin, zin)
      % [X, Z] = SIMULATE(SELF, D, X0, P, Z0)
      %% Allow struct for parameters differential states, algebraic states;
      %% convert to column vectors
      if (isstruct(pin))
	fn = self.model.p;
	for i = 1:numel(fn)
	  p(i) = pin.(fn{i});
	end%for
	p = p(:);
      else
	p = pin;
      end%if
      if (isstruct(xin))
	fn = self.model.x;
	for i = 1:numel(fn)
	  x0(i) = xin.(fn{i});
	end%for
	x0 = x0(:);
      else
	x0 = xin;
      end%if
      % Number of experiments
      nsets = size(x0, 2);

      %% Check dimensions
      assert(size(x0, 1)==self.nx)
      assert(size(p, 1)==self.np)
      assert(size(p, 2)==1)

      % z0 defaults to zero
      if nargin<5
        z0 = zeros(self.nz, nsets);
      else
	if (isstruct(zin))
	  fn = self.model.z;
	  for i = 1:numel(fn)
	    z0(i) = zin.(fn{i});
	  end%for
	  z0 = z0(:);
	else
	  z0 = zin;
	end%if
	assert(size(z0, 1)==self.nz)
	assert(size(z0, 2)==nsets)
      end

      % Solution trajectories
      nt = numel(self.tout);
      x = zeros(self.nx, nt, nsets);
      z = zeros(self.nz, nt, nsets);

      % Simulate for each set of experiments
      for i = 1:nsets
        % Simulate the trajectory for the data set
        sol = self.simulator('x0', x0(:,i), 'p', [p;d], 'z0', z0(:,i));
        x(:,:,i) = full(sol.xf);
        z(:,:,i) = full(sol.zf);
        % Evaluate the measurement function
        sol = self.stagefun('t', self.tout, 'x', sol.xf, 'z', sol.zf,...
                            'p', p, 'd', d);
      end
    end

    function collocation(self)
      % COLLOCATION(SELF) Collocation equations needed for NLP transcription

      % Get collocation coefficients
      [self.tau_root, C, D] = paresto.coll_coeff(self.ord);

      % Parameter vector
      p = casadi.MX.sym('p', self.np);

      % Measurements/outputs
      d = casadi.MX.sym('d', self.nd);

      % Time offset and interval length
      t = casadi.MX.sym('t');
      h = casadi.MX.sym('h');

      % Initial state
      x0 = casadi.MX.sym('x_0', self.nx);

      % State and algebraic variable at collocation points
      x = {};
      z = {};
      xz = {};
      for j=1:self.ord
          x{j} = casadi.MX.sym(['x_' num2str(j)], self.nx);
          xz{end+1} = x{j};
          z{j} = casadi.MX.sym(['z_' num2str(j)], self.nz);
          xz{end+1} = z{j};
      end

      % Expression for the end state
      xf = D(1)*x0;
      for j=1:self.ord
        xf = xf + D(j+1)*x{j};
      end


      % Collocation equations
      g = {};
      for j=1:self.ord
         % Expression for the state derivative
         xp = C(1,j+1)*x0;
         for r=1:self.ord
             xp = xp + C(r+1,j+1)*x{r};
         end
         % Evaluate the DAE right-hand-side
         fj = self.daefun('t', t+h*self.tau_root(j), 'p', p,...
                          'x', x{j}, 'z', z{j}, 'd', d);
         % Collect equality constraints
         g{end+1} = h*fj.ode - xp;
         g{end+1} = fj.alg;
      end

      % Concatenate variables, equations
      xz = vertcat(xz{:});
      g = vertcat(g{:});
      switch self.transcription
        case 'simultaneous'
          % Encapsulate equations in a function object
          self.dynfun = casadi.Function('dynfun', ...
                {x0, [t;h;p;d], xz}, {xf, g}, {'x0', 'p', 'xz'}, {'xf', 'g'});
        case 'shooting'
          % Initial algebraic state
          z0 = casadi.MX.sym('z0', self.nz);
          zf = z{self.ord};
          xzf = [xf;zf];
          % Rootfinding problem
          rfp = struct('x', xz, 'p', [x0;t;h;p;d], 'g', g);
          % Rootfinding solver
	  % rf = casadi.rootfinder('rf', 'newton', rfp);
	  % Suggestion to better handle bad initial guess; Joris Gillis, 9/4/2020
          rf = casadi.rootfinder('rf', 'newton', rfp,struct('error_on_fail',false,'line_search',false));
	  % Function that evaluates state at end time
          xzf_fun = casadi.Function('xzf_fun', {xz}, {xzf}, {'xz'}, {'xzf'}, struct('allow_free', true));
          % Solution to the rootfinding problem
          rfsol = rf('x0', [repmat([x0;z0], self.ord, 1)], 'p', [x0;t;h;p;d]);
          xzfsol = xzf_fun(rfsol.x);
          % Wrap rootfinder in an object with integrator syntax
          self.dynfun = casadi.Function('dynfun', {x0, z0, [t;h;p;d]}, {xzfsol}, ...
                      {'x0', 'z0', 'p'}, {'xzf'});
        case 'sundials'
          if self.nz==0
            plugin = 'cvodes';
          else
            plugin = 'idas';
          end%if
          tau = casadi.MX.sym('tau');
          x = casadi.MX.sym('x', self.daefun.sparsity_in('x'));
          z = casadi.MX.sym('z', self.daefun.sparsity_in('z'));
          fj = self.daefun('t', t + tau*h, 'p', p,...
                           'x', x, 'z', z, 'd', d);
          dae = struct('x', x, 'z', z, 't', tau, 'p', [t; h; p; d], ...
                       'ode', fj.ode*h, 'alg', fj.alg);
          self.dynfun = casadi.integrator('dynfun', plugin, dae, ...
                                           self.nlp_solver_integrator_options);
        otherwise
          error(['No such transcription: ' self.transcription]);
      end
    end

    function transcribe(self)
      % TRANSCRIBE(SELF) NLP transcription

      % Start with an empty NLP
      w={};
      lsq = {};
      g={};

      % w with internal variables eliminated
      w_elim = {};

      % x, z, y, d trajectories
      x = {};
      z = {};
      d = {};

      % Which time points to include in least square term
      lsq_ind = zeros(self.N+1, 1);
      if isfield(self.model, 'lsq_ind')
        % Include only specified
        lsq_ind(self.model.lsq_ind) = 1;
      else
        % Include all
        lsq_ind(:) = 1;
	self.model.lsq_ind = lsq_ind;
      end

      % Sensitivity parameters
      sens = {};

      % Free parameter
      p = casadi.MX.sym('p', self.np);
      w{end+1} = p;
      w_elim{end+1} = p;

      % Perturbation in p
      sens{end+1} = p;

      % Time points and interval lengths
      t = self.tout;
      h = t(2:end)-t(1:end-1);

      % For all experiments
      for e=1:self.nsets
        % Name suffix (multiple experiments only)
        if self.nsets>1
          s = ['_' num2str(e)];
        else
          s = '';
        end

        % Time points
        tk = t(1);

        % Initial conditions
        xk = casadi.MX.sym(['x0' s], self.nx);
        w{end+1} = xk;
        x{end+1} = xk;
        w_elim{end+1} = xk;

        % Initial algebraic variable
        zk = casadi.MX.sym(['z0' s], self.nz);
        w{end+1} = zk;
        z{end+1} = zk;
        w_elim{end+1} = zk;

        % Perturbation in initial conditions
        sens{end+1} = xk;

        % Measurements at initial time
        dk = casadi.MX.sym(['d0' s], self.nd);
        d{end+1} = dk;
        sf = self.stagefun('t', tk, 'x', xk, 'z', zk, 'p', p, 'd', dk);
        g{end+1} = sf.alg;
        
        if lsq_ind(1)
          lsq{end+1} = sf.lsq;
        end

        % Loop over measurements
        for k=1:self.N
          % Interval length
          hk = h(k);

          switch self.transcription
            case 'simultaneous'
              % State and algebraic variables at collocation points
              xzc = casadi.MX.sym(['xz_' num2str(k) s], (self.nx+self.nz)*self.ord);
              w{end+1} = xzc;
              xzk = [xk;zk];
              % w_elim{end+1} = repmat(xzk, self.ord, 1);
              w_elim{end+1} = repmat(xzk, self.ord, 1);
              % Evaluate collocation equations
              Fk = self.dynfun('x0', xk, 'xz', xzc, 'p', [tk;hk;p;dk]);
              g{end+1} = Fk.g;
              
            case 'shooting'
              % Call ODE/DAE integrator
              Fk = self.dynfun('x0', xk, 'z0', zk, 'p', [tk;hk;p;dk]);
            case 'sundials'
              Fk = self.dynfun('x0', xk, 'z0', zk, 'p', [tk;hk;p;dk]);
            otherwise
              error(['No such transcription: ' self.transcription]);
          end

          % New NLP variable for state at end of interval
          xk = casadi.MX.sym(['x' num2str(k) s], self.nx);
          % New algebraic variable at the end of the interval
          zk = casadi.MX.sym(['z' num2str(k) s], self.nz);
          % Measurements
          dk = casadi.MX.sym(['d' num2str(k) s], self.nd);
          tk = t(k+1);
          xzk = [xk;zk];

          sf = self.stagefun('t', tk, 'x', xk, 'z', zk, 'p', p, 'd', dk);
          w{end+1} = xk;
          x{end+1} = xk;
          w_elim{end+1} = xk;

          % Enforce continuity
          switch self.transcription
            case 'simultaneous'
              g{end+1} = xk-Fk.xf;
              g{end+1} = sf.alg;
            case 'shooting'
              g{end+1} = xzk-Fk.xzf;
            case 'sundials'
              g{end+1} = xk-Fk.xf;
              g{end+1} = sf.alg;
            otherwise
              error(['No such transcription: ' self.transcription]);
          end

          w{end+1} = zk;
          z{end+1} = zk;
          w_elim{end+1} = zk;
          d{end+1} = dk;


          if lsq_ind(k+1)
            lsq{end+1} = sf.lsq;
          end
        end
      end

      % Concatenate vectors
      g = vertcat(g{:});
      w = vertcat(w{:});
      sens = vertcat(sens{:});
      lsq = vertcat(lsq{:});
      w_elim = vertcat(w_elim{:});

      % Number of NLP decision variables
      self.nw = numel(w);

      % Sum-of-squares objective
      J = sumsqr(lsq);

      % Trajectories
      x = horzcat(x{:});
      z = horzcat(z{:});
      d = horzcat(d{:});

      % Mapping to and from NLP decision vector
      self.to_w = casadi.Function('to_w', {x, z, p}, {w_elim}, {'x', 'z', 'p'}, {'w'});
      self.from_w = casadi.Function('from_w', {w}, {x, z, p}, {'w'}, {'x', 'z', 'p'});

      % Mapping from sensitivitity parameter to w
      to_sens = casadi.Function('to_sens', {w}, {sens}, {'w'}, {'sens'});
      thetaind = to_sens(1:numel(w));
      self.thetaind = full(thetaind);

      % Create NLP solver
      nlp = struct('f', J, 'x', w, 'g', g, 'p', d);
      self.solver = casadi.nlpsol('solver', self.nlp_solver, nlp, self.nlp_solver_options);

      % If requested, also create a presolver
      if ~isempty(self.nlp_presolver)
        self.presolver = casadi.nlpsol('presolver', self.nlp_presolver, nlp, self.nlp_presolver_options);
      end
    end

    function [r,yy,pp] = optimize(self, d, sol, lb, ub, calc_hess)
      % [R,V] = OPTIMIZE(SELF, D, PHI0, LBPHI, UBPHI, CALC_HESS) Optimize
      % Solve parameter estimation problem and perform a sensitivity analysis
      % of the solution.

      % Get parameters
      p0 = self.struct2vec(sol, 'p', 1, 1, [], 1);
      lbp = self.struct2vec(lb, 'p', 1, 1, -inf, 1);
      ubp = self.struct2vec(ub, 'p', 1, 1, inf, 1);

      % Get state trajectories
      x0 = self.struct2vec(sol, 'x', self.N + 1, self.nsets, [], 1);
      lbx = self.struct2vec(lb, 'x', self.N + 1, self.nsets, -inf, 0);
      ubx = self.struct2vec(ub, 'x', self.N + 1, self.nsets, inf, 0);

      % Get algebraic trajectories
      z0 = self.struct2vec(sol, 'z', self.N + 1, self.nsets, [], 1);
      lbz = self.struct2vec(lb, 'z', self.N + 1, self.nsets, -inf, 0);
      ubz = self.struct2vec(ub, 'z', self.N + 1, self.nsets, inf, 0);
      % Translate to initial guess and bound on w
      w0 = self.to_w(x0, z0, p0);
      lbw = self.to_w(lbx, lbz, lbp);
      ubw = self.to_w(ubx, ubz, ubp);
      traj_indices = setdiff(1:length(w0), self.thetaind);
      lbw(traj_indices) = -inf;
      ubw(traj_indices) = inf;
      % Log message with timings
      if (self.print_level > 0)
        msg = @(m) fprintf('paresto.paresto (t=%g ms): %s\n', 1000*toc, m);
	tic;
      else
        msg = @(m) fprintf('');
      end

      % calc_hess true by default
      if nargin<6
        calc_hess = true;
      end

      % Check dimensions
      assert(size(d, 1)==self.nd);
      assert(size(d, 2)==self.N+1);
      assert(size(d, 3)==self.nsets);

      % Flatten d
      d = reshape(d, self.nd, (self.N+1)*self.nsets);

      % Solution guess
      sol = struct('x', w0, 'lam_x', 0, 'lam_g', 0);

      % Presolve the NLP
      if ~isempty(self.nlp_presolver)
        msg('Presolving NLP');
        sol = self.presolver('x0', sol.x, 'lam_x0', sol.lam_x, 'lam_g0', sol.lam_g,...
                             'lbx', lbw, 'ubx', ubw, 'lbg', 0, 'ubg', 0, 'p', d);
      end

      % Solve the NLP
      msg('Solving NLP'); 
      x0 = sol.x;
      lam_g0 = sol.lam_g;
      lam_x0 = sol.lam_x;
      %disp(self.solver.get_function('nlp_g'))
      sol = self.solver('x0', x0, 'lam_x0', lam_x0, 'lam_g0', lam_g0,...
                        'lbx', lbw, 'ubx', ubw, 'lbg', 0, 'ubg', 0, 'p', d);
      % disp(sol.x)
      nlp_g_function = self.solver.get_function('nlp_g');
      % disp(class(nlp_g_function));
      % disp(d);
      result = nlp_g_function(sol.x,d);
      % disp(result)
      % Return structure
      r = struct;

      % Get the estimated parameters
      w_opt = full(sol.x);
      r.thetavec = w_opt(self.thetaind);

      % Fix the parameters and resolve the NLP
      msg('Resolving NLP with fixed parameters');
      lbw(self.thetaind) = r.thetavec;
      ubw(self.thetaind) = r.thetavec;


      sol = self.solver('x0', sol.x, 'lam_x0', sol.lam_x, 'lam_g0', sol.lam_g,...
                   'lbx', lbw, 'ubx', ubw, 'lbg', 0, 'ubg', 0, 'p', d);

      % Optimal cost
      r.f = full(sol.f);

      % Get solution trajectories
      [x, z, p] = self.from_w(sol.x);
      
      r.x = reshape(full(x), self.nx, self.N + 1, self.nsets);
      r.z = reshape(full(z), self.nz, self.N + 1, self.nsets);
      r.p = full(p);

      % Split up p by variable name
      pp = self.parfun('p', p);
      fn = fieldnames(pp);
      for i=1:numel(fn)
        pp.(fn{i}) = full(pp.(fn{i}));
      end

      % Split up trajectories by variable name
      yy = self.outfun('t', repmat(self.tout, 1, self.nsets), 'x', x, ...
                       'z', z, 'p', p, 'd', d);
      fn = fieldnames(yy);
      for i=1:numel(fn)
        yy.(fn{i}) = reshape(full(yy.(fn{i})), 1, self.N+1, self.nsets);
      end

      % Get multiplier trajectories
      [lam_x, lam_z, lam_p] = self.from_w(sol.lam_x);
      r.lam_x = reshape(full(lam_x), self.nx, self.N + 1, self.nsets);
      r.lam_z = reshape(full(lam_z), self.nz, self.N + 1, self.nsets);
      r.lam_p = full(lam_p);

      % Sensitivities w.r.t. to measurements
      r.lam_d = reshape(full(sol.lam_p), self.nd, self.N+1, self.nsets);

      % Parametric sensitivities
      lam_w = full(sol.lam_x);
      r.df_dtheta = -lam_w(self.thetaind);

      % Forward seeds
      n_est = numel(self.thetaind);
      seed = zeros(self.nw, n_est);
      for i=1:n_est
        seed(self.thetaind(i), i) = 1;
      end

      % Forward sensitivity analysis
      if calc_hess
        msg('Sensitivity analysis');
  %% Ensure lam_g is not exactly zero
	sol.lam_g = full(sol.lam_g);
	sol.lam_g(sol.lam_g == 0) = 1e-300;
	sol.lam_g = casadi.DM(sol.lam_g);

  sol.lam_x = full(sol.lam_x);
  %% Ensure lam_x(~self.thetaind) is exactly zero, 
  %% i.e., trajectory variables x1, x2,... and z0, z1,... are unconstrained
  sol.lam_x(traj_indices) = 0;
  %% Ensure lam_x(self.thetaind) is not exactly zero,
  %% i.e., enforce equality constraints 
  %% on estimated parameters and initial conditions
  %% even if they don't appear in objective function
  sol.lam_x(sol.lam_x(self.thetaind)==0) = 1e-300;
  sol.lam_x = casadi.DM(sol.lam_x);



	fsol = self.fsolver('x0', sol.x, 'lam_x0', sol.lam_x, 'lam_g0', sol.lam_g,...
			    'lbx', lbw, 'ubx', ubw, 'lbg', 0, 'ubg', 0, 'p', d,...
			    'out_x', sol.x, 'out_lam_x', sol.lam_x, 'out_lam_g', sol.lam_g,...
			    'out_lam_p', sol.lam_p, 'out_f', sol.f, 'out_g', sol.g,... 
			    'fwd_lbx', seed, 'fwd_ubx', seed);
        sens = -full(fsol.fwd_lam_x);
        r.d2f_dtheta2 = sens(self.thetaind,:);

        % Get forward derivatives w.r.t. theta
        [dx_dtheta, dz_dtheta, dp_dtheta] = self.from_w(fsol.fwd_x);
        r.dx_dtheta = reshape(full(dx_dtheta), self.nx, self.N + 1, self.nsets, n_est);
        r.dz_dtheta = reshape(full(dz_dtheta), self.nz, self.N + 1, self.nsets, n_est);
        r.dp_dtheta = full(dp_dtheta);

        % Fields in theta
        r.thetafields = self.model.p;
	for i = 1:self.np
	  field = self.model.p{i};
	  lbnames.(field) = lb.(field);
	  ubnames.(field) = ub.(field);
	end
        for e=1:self.nsets
          % Name suffix (multiple experiments only)
          if self.nsets>1
            s = ['_' num2str(e)];
          else
            s = '';
          end
          % Add initial condition fields; create names for IC bounds
          for i=1:self.nx
	    % Initial condition field name
	    ICname = [self.model.x{i} '0' s];
	    r.thetafields{end+1} = ICname;
	    if (isfield (lb, self.model.x{i}))
	      if self.nsets > 1
		lbnames.(ICname) = lb.(self.model.x{i})(:,:,e);
		ubnames.(ICname) = ub.(self.model.x{i})(:,:,e);
	      else
		lbnames.(ICname) = lb.(self.model.x{i})(:,1);
		ubnames.(ICname) = ub.(self.model.x{i})(:,1);
              end
	    else
	      lbnames.(ICname) = -inf;
	      ubnames.(ICname) = inf;
	    end
	  end
        end
        % Split up dx_dtheta, dz_dtheta by name
        for j=1:numel(r.thetafields)
          for i=1:self.nx
            r.(['d' self.model.x{i} '_d' r.thetafields{j}]) = ...
              reshape(r.dx_dtheta(i, :, :, j), self.N + 1, self.nsets);
          end
          for i=1:self.nz
            r.(['d' self.model.z{i} '_d' r.thetafields{j}]) = ...
              reshape(r.dz_dtheta(i, :, :, j), self.N + 1, self.nsets);
          end
        end
      end
      % Return all estimated parameters by name
      for j = 1: numel(r.thetafields)
	r.par.(r.thetafields{j}) = r.thetavec(j);
      end
      % Store index of estimated parameters not having equality constraints
      r.conf_ind = find(cell2mat(struct2cell(lbnames)) < cell2mat(struct2cell(ubnames)));
      % Return by name estimated parameters not having equality constraints
      for i = 1: numel(r.conf_ind)
	r.theta.(r.thetafields{r.conf_ind(i)}) = r.thetavec(r.conf_ind(i));
      end
      %% Store best guess of number of data points in each dataset:
      %% number time points x length of measurement vector
      r.n_data = numel(self.model.lsq_ind)*numel(self.model.d);
      % Done
      msg('Optimization complete');
    end

    function conf = confidence(self, r, alpha)
%    function theta_conf = confidence(self, r, conf_ind, alpha)
      % THETA_CONF = CONFIDENCE(R, CONF_IND, ALPHA): Calculate confidence intervals

      % conf_ind, index of parameters not having equality constraints,
      % defined in optimize function
      conf_ind = r.conf_ind;

      % alpha defaults to 0.95
      if nargin<3
        alpha = 0.95;
      end

      % Number of parameters being estimated
      n_est = numel(conf_ind);

      % Quick return if n_est = 0
      if n_est==0
	conf.H =[];
	conf.diag_inv_H = [];
	conf.bbox = [];
	conf.mbox = [];
        return
      end

      % Get the subset of the reduced Hessian being estimated
      H = r.d2f_dtheta2(conf_ind, conf_ind);
      % Ensure symmetry
      if norm(H-H.', 'inf')>1e-6*norm(H, 'inf')
        warning('Reduced Hessian appears nonsymmetric');
      end
      H = 0.5*(H.' + H);
      [v,e] = eig(H);
      e = diag(e);
      %%% Inspect the eigenvalues, recursive call if non-positive eigenvalues
      %% if (any(e<1e-10))
      %%   theta_conf = inf(n_est, 1);
      %%   i = find(e>=1e-10);
      %%   % Call recursively
      %%   theta_conf(i) = self.confidence(r, conf_ind(i), alpha);
      %%   return;
      %% end
      rH = sum ( e >= 1e-10 );
      if( rH < n_est )
	%% treat small or negative eigenvalues here
	[e, perm] = sort(e, 'descend');
	v = v(:, perm);
	v1 = v(:,1:rH);
	v2 = v(:,rH+1:end);
	diag_inv_H = diag( v1*diag( 1./e(1:rH) )*v1' );
	%% insert Inf for entries diag(inv(H)) where inverting suspect eigenvalue
	ind = (diag(v2*v2') ~= 0);
	diag_inv_H(ind) = inf;
      else
	diag_inv_H = diag(v*diag(1./e)*v');
      end
      % Total number of data points; either user has provided it, or 
      % make educated guess
      if (isfield(self.model, 'ndata'))
        n_data = self.ndata;
      else
        n_data = self.nsets*r.n_data;
      end
      % Calculate Fstat, bounding box and marginal box
      try
        Fstat = finv(alpha, n_est, n_data-n_est);
        Fstatm = finv(alpha, 1, n_data-n_est);
      catch ME
        try
          % Try to load finv from a file named e.g. finv95.mat for alpha=0.95
          finv_mat = ['finv' num2str(alpha*100)];
          S = load(finv_mat);
          M = S.(finv_mat);
          % Get the value as the corresponding matrix entry
          Fstat = full(M(n_est, n_data-n_est));
	  Fstatm = full(M(1, n_data-n_est));
          % Make sure it's not zero (e.g. missing entry in sparse matrix)
          assert(Fstat~=0, 'Entry not available');
        catch ME2
          % Informative error message
          warning(ME.message);
          warning(ME2.message);
          error('paresto:fstat', ['finv is not available. ',...
              'Try to manually provide model.Fstat=finv(%g,%d,%g) or a ',...
              'file named ' finv_file '.'], alpha, n_est, n_data-n_est);
        end
      end

      % Return Hessian, diagonal of inverse, and confidence intervals (bounding and marginal)
      conf.H =H;
      conf.diag_inv_H = diag_inv_H;
      %inv_hess = inv(H);
      theta_conf = sqrt(2*n_est/(n_data-n_est)*Fstat*r.f*diag_inv_H);
      theta_marg = sqrt(Fstatm/(Fstat*n_est))*theta_conf;
      % Return confidence intervals by name
      for i = 1: n_est
	fn = r.thetafields{r.conf_ind(i)};
	conf.bbox.(fn) = theta_conf(i);
	conf.mbox.(fn) = theta_marg(i);
      end
    end

    function [a,v] = str2sym(self, fname, v)
      % [A,V] = STR2SYM(SELF,FNAME,V) Create CasADi symbols
      assert(isfield(self.model, fname))
      s = self.model.(fname);
      assert(iscell(s));
      n = numel(s);
      % Quick return
      if n==0
        a = casadi.SX(0,1);
        return
      end
      % Create symbols
      a = cell(1,n);
      for i=1:n
        % Component name
        si = s{i};

        % Get dimension
        if isfield(self.model, 'dim') && isfield(self.model.dim, si)
          d = self.model.dim.(si);
        else
          d = 1; % scalar by default
        end

        % Create symbol
        a{i} = casadi.SX.sym(si, d);
        % Make sure that it does not already exist in v
        assert(~isfield(v, si), ['Duplicate expression: ' si]);
        v.(si) = a{i};
      end
      a = vertcat(a{:});
    end

    function a = fun2sym(self, fname, t, y, p)
      % A = FUN2SYM(SELF,FNAME,T,Y,P) Create symbols from function handle
      % Quick return if not provided
      if ~isfield(self.model, fname)
        r = struct;
        a = casadi.SX(0, 1);
        return;
      end
      % Create struct
      try
        % User provided function
        a = self.model.(fname)(t, y, p);
        assert(iscell(a), 'Expected cell output')
        a = vertcat(a{:});
      catch ME
        % More informative error message
        warning(ME.message);
        error(['paresto:' fname],...
              'Failure evaluating model.%s(t, y, p).\n%s\n%s', fname, ...
              ['y has fields ', strjoin(fieldnames(y), ',')], ...
              ['p has fields ', strjoin(fieldnames(p), ',')]);
      end
    end

    function v = struct2vec(self, s, fname, nrhs, nsets, def, copy)
      % V = STRUCT2VEC(SELF,S,FNAME,NRHS) Get vector from structure

      % Default argument
      assert(isempty(def) || numel(def)==1);

      assert(isfield(self.model, fname))
      ff = self.model.(fname);
      assert(iscell(ff));
      n = numel(ff);
      % Quick return
      if n==0
        v = [];
        return
      end
      % Create symbols
      v = cell(1,n);
      for i=1:n
        % Component name
        f = ff{i};
        % Get dimension
        if isfield(self.model, 'dim') && isfield(self.model.dim, f)
          d = self.model.dim.(f);
        else
          d = 1; % scalar by default
        end
        % Get value
        if isfield(s, f)
          vi = s.(f);
        elseif isempty(def)
          error('%s has not been provided', f);
        else
          vi = def; % default value
        end
        % Make sure at most 3 dimensions
        assert(ndims(vi)<=3);
        % Correct number of elements if needed
        if size(vi, 1)~=d
          assert(mod(d, size(vi, 1))==0);
          vi = repmat(vi, d/size(vi, 1), 1);
	end
	% Correct number of right-hand sides if needed
	if size(vi, 2)~=nrhs
          assert(mod(nrhs, size(vi, 2))==0);
	  %% fill right-hand sides with copy of given vi
	  fill = 1+size(vi,2):nrhs;
	  vi = repmat(vi, 1, nrhs/size(vi, 2));
          if (~copy)
	    %% -or- overwrite right-hand sides with default
	    vi(:,fill,:) = def;
	  end
        end
        % Correct number of data sets if needed
        if size(vi, 3)~=nsets
          assert(mod(nsets, size(vi, 3))==0);
          vi = repmat(vi, 1, 1, nsets/size(vi, 3));
        end
        % Flatten second and third dimensions, if needed
        if size(vi, 3)~=1
          vi = reshape(vi, size(vi, 1), []);
        end
        % Save to list
        v{i} = vi;
      end
      % Concatenate
      v = vertcat(v{:});
    end
  end

  methods(Static = true)
    function [tau_root, C, D, B] = coll_coeff(ord)
      % [TAU_ROOT, C, D, B] = COLL_COEFF(ORD) Get collocation coefficients

      % Gauss-Radau points
      switch ord
        case 1
          % tau_root = [0.5];
          tau_root = [1.];
        case 2
          % tau_root = [0.5-sqrt(3)/6; 0.5+sqrt(3)/6];
          tau_root = [1/3; 1];
        case 3
        %   tau_root = [0.5-sqrt(15)/10; 0.5; 0.5+sqrt(15)/10];
          tau_root = [2/5 - sqrt(6)/10; 2/5 + sqrt(6)/10; 1];
        otherwise
          error('order not supported');
       end

       % Degree of interpolating polynomial
       d = numel(tau_root);

       % Augment with zero
       t = [0; tau_root];

       % Coefficients of the collocation equation
       C = zeros(ord+1,ord+1);

       % Coefficients of the continuity equation
       D = zeros(ord+1, 1);

       % Coefficients of the quadrature function
       B = zeros(ord+1, 1);

       % Construct polynomial basis
       for j=1:ord+1
         % Construct Lagrange polynomials to get the polynomial basis at the collocation point
         coeff = 1;
         for r=1:ord+1
           if r ~= j
             coeff = conv(coeff, [1, -t(r)]);
             coeff = coeff / (t(j)-t(r));
           end
         end
         % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
         D(j) = polyval(coeff, 1.0);

         % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
         pder = polyder(coeff);
         for r=1:ord+1
           C(j,r) = polyval(pder, t(r));
         end

         % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
         pint = polyint(coeff);
         B(j) = polyval(pint, 1.0);
       end
    end

    function generate_finv(alpha, m_range, n_range)
      % GENERATE_FINV(D) Generate a file with finv values
      % For use when finv is not available
      [M,N] = meshgrid(m_range, n_range);
      Z = finv(alpha,M,N);
      % Place in a large enough sparse matrix
      M = sparse(max(m_range), max(n_range));
      M(m_range, n_range) = Z;
      % Name of the variable and file
      finv_mat = ['finv' num2str(alpha*100)];
      % Save to a file
      S = struct(finv_mat, M);
      save('-mat', [finv_mat '.mat'], '-struct', 'S', finv_mat);
    end

    function [x, y, major, minor, bbox] = ellipse (amat, level, n, shift)
    % [x, y, major, minor, bbox] = ELLIPSE(amat, level, n, shift)
    %
    % Given a 2x2 matrix, generate ellipse data for plotting.  The
    % arguments N and SHIFT are optional.  If N is an empty matrix, a
    % default value of 100 is used.
    % James B. Rawlings and John W. Eaton, 2001

      % Check number of arguments
      narginchk(2, 4);

      % Set default arguments
      if (nargin < 3)
        n = 100;
        if (nargin < 4)
          shift = [0, 0];
        end
      end
      if (isempty(n))
        n = 100;
      end

      ss = size(shift);
      if (any(ss ~= [1, 2]))
        if (ss == [2, 1])
          shift = shift';
        else
          error('shift must be a 2-element row vector');
        end
      end

      [v, l] = eig(amat / level);
      dl = diag(l);
      if (any(imag (dl)) || any(dl <= 0))
        error('ellipse: amat must be positive definite');
      end

      % Generate contour data.
      a = 1 / sqrt(l(1,1));
      b = 1 / sqrt(l(2,2));
      t = linspace(0, 2*pi, n)';
      xt = a * cos(t);
      yt = b * sin(t);

      % Rotate the contours.
      ra = atan2(v(2,1), v(1,1));
      cos_ra = cos(ra);
      sin_ra = sin(ra);
      x = xt * cos_ra - yt * sin_ra + shift(1);
      y = xt * sin_ra + yt * cos_ra + shift(2);

      % Endpoints of the major and minor axes.
      minor = (v * diag([a, b]))';
      major = minor;
      major(2,:) = -major(1,:);
      minor(1,:) = -minor(2,:);
      t = [1; 1] * shift;
      major = major + t;
      minor = minor + t;

      % Bounding box for the ellipse using magic formula.
      ainv = inv(amat);
      xbox = sqrt(level * ainv(1,1));
      ybox = sqrt(level * ainv(2,2));
      bbox = [xbox ybox; xbox -ybox; -xbox -ybox; -xbox ybox; xbox ybox];
      t = [1; 1; 1; 1; 1] * shift;
      bbox = bbox + t;
    end
  end
end
