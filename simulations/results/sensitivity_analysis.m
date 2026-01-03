
% Sensitivity Analysis using Monte Carlo
% Variation of geometric and material parameters

n_samples = 100;  % Number of Monte Carlo samples

% Storage for results
sensitivity_results = struct();
sensitivity_results.samples = n_samples;
sensitivity_results.thrust = zeros(n_samples, 1);
sensitivity_results.energy = zeros(n_samples, 1);
sensitivity_results.max_temp = zeros(n_samples, 1);

% Nominal parameter values
nominal = struct();
nominal.r_nozzle_inner = 5e-05;
nominal.gap_height = 0.0002;
nominal.r_toroid_minor = 0.0003;
nominal.pitch_coil = 0.0002;
nominal.Mo_conductivity = 18900000.0;
nominal.AlN_d33 = 5.53e-12;
nominal.Sapphire_n = 1.76;

% Monte Carlo loop
for i = 1:n_samples
    % Perturb parameters (normal distribution)
    model.param.set('r_nozzle_inner', nominal.r_nozzle_inner * (1 + 0.1 * randn()));
    model.param.set('gap_height', nominal.gap_height * (1 + 0.2 * randn()));
    model.param.set('r_toroid_minor', nominal.r_toroid_minor * (1 + 0.15 * randn()));
    model.param.set('pitch_coil', nominal.pitch_coil * (1 + 0.1 * randn()));
    model.param.set('Mo_conductivity', nominal.Mo_conductivity * (1 + 0.05 * randn()));
    model.param.set('AlN_d33', nominal.AlN_d33 * (1 + 0.1 * randn()));
    model.param.set('Sapphire_n', nominal.Sapphire_n * (1 + 0.02 * randn()));

    % Run simulation
    model.study('std2').run();  % Frequency domain
    
    % Extract metrics
    model.result.numerical('int1').setIndex('expr', 'emw.Poavz', 0);
    model.result.numerical('int1').setResult();
    sensitivity_results.thrust(i) = model.result.numerical('int1').getReal();
    
    % Energy
    sensitivity_results.energy(i) = mphglobal(model, 'emw.intWe', 'dataset', 'dset2');
    
    % Max temperature
    sensitivity_results.max_temp(i) = mphmax(model, 'T', 'dataset', 'dset2');
    
    if mod(i, 10) == 0
        fprintf('Completed %d/%d samples\n', i, n_samples);
    end
end

% Statistical analysis
fprintf('\n=== SENSITIVITY ANALYSIS RESULTS ===\n');
fprintf('Thrust: mean = %e, std = %e (CV = %.2f%%)\n', ...
    mean(sensitivity_results.thrust), std(sensitivity_results.thrust), ...
    100*std(sensitivity_results.thrust)/mean(sensitivity_results.thrust));

% Save results
save('results/sensitivity_analysis.mat', 'sensitivity_results');
