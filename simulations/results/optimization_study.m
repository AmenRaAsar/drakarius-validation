
% Optimization Study: Maximize Thrust
% Optimize: plasma density, actuation voltage, coil current
% Constraint: Temperature < 300K, Power < 10kW

% Define optimization variables
opt_vars = {'n_plasma', 'V_aln', 'I_coil'};

% Bounds
lower_bounds = [1e16, 50, 1];     % [m^-3, V, A]
upper_bounds = [1e18, 200, 20];   % [m^-3, V, A]

% Initial guess (nominal values)
x0 = [1e17, 100, 10];

% Objective function (negative thrust for minimization)
objective = @(x) -compute_thrust(model, x);

% Constraint functions
nonlcon = @(x) optimization_constraints(model, x);

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6);

% Run optimization
[x_opt, thrust_opt] = fmincon(objective, x0, [], [], [], [], ...
    lower_bounds, upper_bounds, nonlcon, options);

fprintf('\n=== OPTIMIZATION RESULTS ===\n');
fprintf('Optimal plasma density: %e m^-3\n', x_opt(1));
fprintf('Optimal actuation voltage: %.1f V\n', x_opt(2));
fprintf('Optimal coil current: %.1f A\n', x_opt(3));
fprintf('Maximum thrust: %e N\n', -thrust_opt);

% Save optimal configuration
optimal_config = struct();
optimal_config.n_plasma = x_opt(1);
optimal_config.V_aln = x_opt(2);
optimal_config.I_coil = x_opt(3);
optimal_config.thrust = -thrust_opt;

save('results/optimal_configuration.mat', 'optimal_config');

function thrust = compute_thrust(model, x)
    % Set parameters
    model.param.set('n_plasma', sprintf('%e[1/m^3]', x(1)));
    model.param.set('V_aln', sprintf('%g[V]', x(2)));
    model.param.set('I_coil', sprintf('%g[A]', x(3)));
    
    % Run simulation
    model.study('std2').run();
    
    % Compute thrust
    model.result.numerical('int1').setIndex('expr', 'emw.Poavz', 0);
    model.result.numerical('int1').setResult();
    thrust = model.result.numerical('int1').getReal();
end

function [c, ceq] = optimization_constraints(model, x)
    % Inequality constraints: c <= 0
    % Equality constraints: ceq = 0
    
    % Set parameters
    model.param.set('n_plasma', sprintf('%e[1/m^3]', x(1)));
    model.param.set('V_aln', sprintf('%g[V]', x(2)));
    model.param.set('I_coil', sprintf('%g[A]', x(3)));
    
    % Run simulation
    model.study('std2').run();
    
    % Temperature constraint: max(T) < 300K
    max_temp = mphmax(model, 'T', 'dataset', 'dset2');
    c(1) = max_temp - 300;
    
    % Power constraint: total power < 10kW
    total_power = mphglobal(model, 'emw.Qh', 'dataset', 'dset2');
    c(2) = total_power - 10e3;
    
    % No equality constraints
    ceq = [];
end
