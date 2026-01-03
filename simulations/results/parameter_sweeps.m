% All Parameter Sweeps for Drakarius Validation


% Parameter Sweep: Plasma density sweep to find optimal operating point
% Parameter: n_plasma

% Add to COMSOL model after creating all studies
model.param.set('n_plasma', 'sweep_var', 'Plasma density sweep to find optimal operating point');

% Create parametric sweep study
model.study.create('param_sweep_plasma_density');
model.study('param_sweep_plasma_density').create('param', 'Parametric');
model.study('param_sweep_plasma_density').feature('param').set('pname', {'n_plasma'});
model.study('param_sweep_plasma_density').feature('param').set('plistarr', {'sweep_values'});
model.study('param_sweep_plasma_density').feature('param').set('punit', {'1/m^3'});

% Sweep values
sweep_values = [1.000000e+16, 1.274275e+16, 1.623777e+16, 2.069138e+16, 2.636651e+16, 3.359818e+16, 4.281332e+16, 5.455595e+16, 6.951928e+16, 8.858668e+16, 1.128838e+17, 1.438450e+17, 1.832981e+17, 2.335721e+17, 2.976351e+17, 3.792690e+17, 4.832930e+17, 6.158482e+17, 7.847600e+17, 1.000000e+18];

% Attach to frequency domain study
model.study('param_sweep_plasma_density').create('freq', 'Frequency');
model.study('param_sweep_plasma_density').feature('freq').set('plist', 'f_plasma');

% Run sweep
model.study('param_sweep_plasma_density').run();

% Extract results
thrust_vs_plasma_density = [];
for i = 1:length(sweep_values)
    model.result.numerical('int1').setIndex('expr', 'emw.Poavz', 0);
    model.result.numerical('int1').setResult();
    thrust_vs_plasma_density(i) = model.result.numerical('int1').getReal();
end

% Save results
save('results/sweep_plasma_density.mat', 'sweep_values', 'thrust_vs_plasma_density');

fprintf('Parameter sweep complete: plasma_density\n');
fprintf('Optimal n_plasma = %g 1/m^3\n', sweep_values(find(thrust_vs_plasma_density == max(thrust_vs_plasma_density))));

% ======================================================================


% Parameter Sweep: AlN actuation voltage sweep for piezoelectric optimization
% Parameter: V_aln

% Add to COMSOL model after creating all studies
model.param.set('V_aln', 'sweep_var', 'AlN actuation voltage sweep for piezoelectric optimization');

% Create parametric sweep study
model.study.create('param_sweep_actuation_voltage');
model.study('param_sweep_actuation_voltage').create('param', 'Parametric');
model.study('param_sweep_actuation_voltage').feature('param').set('pname', {'V_aln'});
model.study('param_sweep_actuation_voltage').feature('param').set('plistarr', {'sweep_values'});
model.study('param_sweep_actuation_voltage').feature('param').set('punit', {'V'});

% Sweep values
sweep_values = [5.000000e+01, 6.000000e+01, 7.000000e+01, 8.000000e+01, 9.000000e+01, 1.000000e+02, 1.100000e+02, 1.200000e+02, 1.300000e+02, 1.400000e+02, 1.500000e+02, 1.600000e+02, 1.700000e+02, 1.800000e+02, 1.900000e+02, 2.000000e+02];

% Attach to frequency domain study
model.study('param_sweep_actuation_voltage').create('freq', 'Frequency');
model.study('param_sweep_actuation_voltage').feature('freq').set('plist', 'f_plasma');

% Run sweep
model.study('param_sweep_actuation_voltage').run();

% Extract results
thrust_vs_actuation_voltage = [];
for i = 1:length(sweep_values)
    model.result.numerical('int1').setIndex('expr', 'emw.Poavz', 0);
    model.result.numerical('int1').setResult();
    thrust_vs_actuation_voltage(i) = model.result.numerical('int1').getReal();
end

% Save results
save('results/sweep_actuation_voltage.mat', 'sweep_values', 'thrust_vs_actuation_voltage');

fprintf('Parameter sweep complete: actuation_voltage\n');
fprintf('Optimal V_aln = %g V\n', sweep_values(find(thrust_vs_actuation_voltage == max(thrust_vs_actuation_voltage))));

% ======================================================================


% Parameter Sweep: Coil current sweep for magnetic field strength
% Parameter: I_coil

% Add to COMSOL model after creating all studies
model.param.set('I_coil', 'sweep_var', 'Coil current sweep for magnetic field strength');

% Create parametric sweep study
model.study.create('param_sweep_coil_current');
model.study('param_sweep_coil_current').create('param', 'Parametric');
model.study('param_sweep_coil_current').feature('param').set('pname', {'I_coil'});
model.study('param_sweep_coil_current').feature('param').set('plistarr', {'sweep_values'});
model.study('param_sweep_coil_current').feature('param').set('punit', {'A'});

% Sweep values
sweep_values = [1.000000e+00, 2.000000e+00, 3.000000e+00, 4.000000e+00, 5.000000e+00, 6.000000e+00, 7.000000e+00, 8.000000e+00, 9.000000e+00, 1.000000e+01, 1.100000e+01, 1.200000e+01, 1.300000e+01, 1.400000e+01, 1.500000e+01, 1.600000e+01, 1.700000e+01, 1.800000e+01, 1.900000e+01, 2.000000e+01];

% Attach to frequency domain study
model.study('param_sweep_coil_current').create('freq', 'Frequency');
model.study('param_sweep_coil_current').feature('freq').set('plist', 'f_plasma');

% Run sweep
model.study('param_sweep_coil_current').run();

% Extract results
thrust_vs_coil_current = [];
for i = 1:length(sweep_values)
    model.result.numerical('int1').setIndex('expr', 'emw.Poavz', 0);
    model.result.numerical('int1').setResult();
    thrust_vs_coil_current(i) = model.result.numerical('int1').getReal();
end

% Save results
save('results/sweep_coil_current.mat', 'sweep_values', 'thrust_vs_coil_current');

fprintf('Parameter sweep complete: coil_current\n');
fprintf('Optimal I_coil = %g A\n', sweep_values(find(thrust_vs_coil_current == max(thrust_vs_coil_current))));

% ======================================================================


% Parameter Sweep: Temperature sweep from cryogenic to room temperature
% Parameter: T_cryo

% Add to COMSOL model after creating all studies
model.param.set('T_cryo', 'sweep_var', 'Temperature sweep from cryogenic to room temperature');

% Create parametric sweep study
model.study.create('param_sweep_temperature');
model.study('param_sweep_temperature').create('param', 'Parametric');
model.study('param_sweep_temperature').feature('param').set('pname', {'T_cryo'});
model.study('param_sweep_temperature').feature('param').set('plistarr', {'sweep_values'});
model.study('param_sweep_temperature').feature('param').set('punit', {'K'});

% Sweep values
sweep_values = [4.000000e+00, 1.000000e+01, 2.000000e+01, 4.000000e+01, 7.700000e+01, 1.500000e+02, 3.000000e+02];

% Attach to frequency domain study
model.study('param_sweep_temperature').create('freq', 'Frequency');
model.study('param_sweep_temperature').feature('freq').set('plist', 'f_plasma');

% Run sweep
model.study('param_sweep_temperature').run();

% Extract results
thrust_vs_temperature = [];
for i = 1:length(sweep_values)
    model.result.numerical('int1').setIndex('expr', 'emw.Poavz', 0);
    model.result.numerical('int1').setResult();
    thrust_vs_temperature(i) = model.result.numerical('int1').getReal();
end

% Save results
save('results/sweep_temperature.mat', 'sweep_values', 'thrust_vs_temperature');

fprintf('Parameter sweep complete: temperature\n');
fprintf('Optimal T_cryo = %g K\n', sweep_values(find(thrust_vs_temperature == max(thrust_vs_temperature))));

% ======================================================================


% Parameter Sweep: Plasma frequency sweep for resonance identification
% Parameter: f_plasma

% Add to COMSOL model after creating all studies
model.param.set('f_plasma', 'sweep_var', 'Plasma frequency sweep for resonance identification');

% Create parametric sweep study
model.study.create('param_sweep_frequency');
model.study('param_sweep_frequency').create('param', 'Parametric');
model.study('param_sweep_frequency').feature('param').set('pname', {'f_plasma'});
model.study('param_sweep_frequency').feature('param').set('plistarr', {'sweep_values'});
model.study('param_sweep_frequency').feature('param').set('punit', {'Hz'});

% Sweep values
sweep_values = [1.000000e+09, 1.100000e+09, 1.200000e+09, 1.300000e+09, 1.400000e+09, 1.500000e+09, 1.600000e+09, 1.700000e+09, 1.800000e+09, 1.900000e+09, 2.000000e+09, 2.100000e+09, 2.200000e+09, 2.300000e+09, 2.400000e+09, 2.500000e+09, 2.600000e+09, 2.700000e+09, 2.800000e+09, 2.900000e+09, 3.000000e+09, 3.100000e+09, 3.200000e+09, 3.300000e+09, 3.400000e+09, 3.500000e+09, 3.600000e+09, 3.700000e+09, 3.800000e+09, 3.900000e+09, 4.000000e+09, 4.100000e+09, 4.200000e+09, 4.300000e+09, 4.400000e+09, 4.500000e+09, 4.600000e+09, 4.700000e+09, 4.800000e+09, 4.900000e+09, 5.000000e+09];

% Attach to frequency domain study
model.study('param_sweep_frequency').create('freq', 'Frequency');
model.study('param_sweep_frequency').feature('freq').set('plist', 'f_plasma');

% Run sweep
model.study('param_sweep_frequency').run();

% Extract results
thrust_vs_frequency = [];
for i = 1:length(sweep_values)
    model.result.numerical('int1').setIndex('expr', 'emw.Poavz', 0);
    model.result.numerical('int1').setResult();
    thrust_vs_frequency(i) = model.result.numerical('int1').getReal();
end

% Save results
save('results/sweep_frequency.mat', 'sweep_values', 'thrust_vs_frequency');

fprintf('Parameter sweep complete: frequency\n');
fprintf('Optimal f_plasma = %g Hz\n', sweep_values(find(thrust_vs_frequency == max(thrust_vs_frequency))));

% ======================================================================

