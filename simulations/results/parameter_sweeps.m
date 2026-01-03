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
sweep_values = [np.float64(1e+16), np.float64(1.2742749857031322e+16), np.float64(1.6237767391887176e+16), np.float64(2.0691380811147816e+16), np.float64(2.6366508987303664e+16), np.float64(3.359818286283788e+16), np.float64(4.281332398719396e+16), np.float64(5.4555947811685144e+16), np.float64(6.951927961775591e+16), np.float64(8.858667904100795e+16), np.float64(1.128837891684693e+17), np.float64(1.438449888287666e+17), np.float64(1.8329807108324374e+17), np.float64(2.3357214690901213e+17), np.float64(2.976351441631313e+17), np.float64(3.792690190732238e+17), np.float64(4.8329302385717325e+17), np.float64(6.158482110660279e+17), np.float64(7.847599703514623e+17), np.float64(1e+18)];

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
sweep_values = [np.float64(50.0), np.float64(60.0), np.float64(70.0), np.float64(80.0), np.float64(90.0), np.float64(100.0), np.float64(110.0), np.float64(120.0), np.float64(130.0), np.float64(140.0), np.float64(150.0), np.float64(160.0), np.float64(170.0), np.float64(180.0), np.float64(190.0), np.float64(200.0)];

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
sweep_values = [np.float64(1.0), np.float64(2.0), np.float64(3.0), np.float64(4.0), np.float64(5.0), np.float64(6.0), np.float64(7.0), np.float64(8.0), np.float64(9.0), np.float64(10.0), np.float64(11.0), np.float64(12.0), np.float64(13.0), np.float64(14.0), np.float64(15.0), np.float64(16.0), np.float64(17.0), np.float64(18.0), np.float64(19.0), np.float64(20.0)];

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
sweep_values = [np.int64(4), np.int64(10), np.int64(20), np.int64(40), np.int64(77), np.int64(150), np.int64(300)];

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
sweep_values = [np.float64(1000000000.0), np.float64(1100000000.0), np.float64(1200000000.0), np.float64(1300000000.0), np.float64(1400000000.0), np.float64(1500000000.0), np.float64(1600000000.0), np.float64(1700000000.0), np.float64(1800000000.0), np.float64(1900000000.0), np.float64(2000000000.0), np.float64(2100000000.0), np.float64(2200000000.0), np.float64(2300000000.0), np.float64(2400000000.0), np.float64(2500000000.0), np.float64(2600000000.0), np.float64(2700000000.0), np.float64(2800000000.0), np.float64(2900000000.0), np.float64(3000000000.0), np.float64(3100000000.0), np.float64(3200000000.0), np.float64(3300000000.0), np.float64(3400000000.0), np.float64(3500000000.0), np.float64(3600000000.0), np.float64(3700000000.0), np.float64(3800000000.0), np.float64(3900000000.0), np.float64(4000000000.0), np.float64(4100000000.0), np.float64(4200000000.0), np.float64(4300000000.0), np.float64(4400000000.0), np.float64(4500000000.0), np.float64(4600000000.0), np.float64(4700000000.0), np.float64(4800000000.0), np.float64(4900000000.0), np.float64(5000000000.0)];

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

