
% Transient Startup Analysis
% Study system behavior from cold start to steady state

% Extended time range for startup
model.study.create('startup');
model.study('startup').create('time', 'Transient');
model.study('startup').label('Transient Startup Analysis');
model.study('startup').feature('time').set('tlist', 'range(0,1e-9,1e-6)');  % 0 to 1 µs

% Initial conditions (cold, no fields)
model.component('comp1').physics('emw').feature('init1').set('Ex', 0);
model.component('comp1').physics('emw').feature('init1').set('Ey', 0);
model.component('comp1').physics('emw').feature('init1').set('Ez', 0);
model.component('comp1').physics('ht').feature('init1').set('T', '4[K]');

% Run startup simulation
fprintf('Running transient startup analysis...\n');
model.study('startup').run();

% Extract time evolution
times = mphglobal(model, 't', 'dataset', 'dset_startup', 'outersolnum', 'all');
thrust_vs_time = mphglobal(model, 'emw.Poavz', 'dataset', 'dset_startup', 'outersolnum', 'all');
energy_vs_time = mphglobal(model, 'emw.intWe', 'dataset', 'dset_startup', 'outersolnum', 'all');

% Calculate startup time (95% of steady-state)
steady_state_thrust = thrust_vs_time(end);
startup_idx = find(thrust_vs_time >= 0.95*steady_state_thrust, 1);
startup_time = times(startup_idx);

fprintf('\n=== STARTUP ANALYSIS RESULTS ===\n');
fprintf('Startup time (95%% steady-state): %.2e s\n', startup_time);
fprintf('Steady-state thrust: %e N\n', steady_state_thrust);

% Save results
startup_results = struct();
startup_results.times = times;
startup_results.thrust = thrust_vs_time;
startup_results.energy = energy_vs_time;
startup_results.startup_time = startup_time;

save('results/startup_analysis.mat', 'startup_results');
