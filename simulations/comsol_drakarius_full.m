% COMSOL Multiphysics Model for Drakarius Propulsion System (MATLAB)
% Comprehensive simulation of plasma-polariton-piezoelectric propulsion
% 
% Physics Modules:
% - Eigenmode analysis for H2 plasma spiral modes
% - Electromagnetic wave propagation
% - Plasma fluid dynamics
% - Heat transfer with cryogenic cooling (4K)
% - Piezoelectric actuation (AlN)
% - Harmonic studies on odd nodes (1,3,5,7...270°)
% - Casimir effect in whispering gallery modes
% - Polariton formation and avoided crossings
% - Poynting vector directional bias analysis
% 
% Author: Drakarius Framework Validation Team

function model = comsol_drakarius_full()

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('DrakariusFullSystem');
model.modelPath('/home/runner/work/drakarius-validation/drakarius-validation/simulations');
model.label('Drakarius Propulsion System - Full Multiphysics');

%% Geometric Parameters
model.param.set('r_nozzle_inner', '50[um]', 'Mo nozzle inner radius');
model.param.set('r_nozzle_outer', '100[um]', 'Mo nozzle outer radius');
model.param.set('h_nozzle', '500[um]', 'Mo nozzle height');
model.param.set('gap_height', '200[um]', 'Gap to toroid');
model.param.set('r_toroid_major', '2[mm]', 'CTE toroid major radius');
model.param.set('r_toroid_minor', '300[um]', 'CTE toroid minor radius');
model.param.set('h_coil', '5[mm]', 'Height of spiral coils');
model.param.set('r_coil', '1.5[mm]', 'Radius of spiral coils');
model.param.set('pitch_coil', '200[um]', 'Coil pitch');
model.param.set('r_resonator_outer', '3[mm]', 'Sapphire resonator outer radius');
model.param.set('r_resonator_inner', '2.5[mm]', 'Sapphire resonator inner radius');
model.param.set('h_resonator', '4[mm]', 'Resonator height');

%% Physics Parameters
model.param.set('T_cryo', '4[K]', 'Cryogenic temperature');
model.param.set('f_plasma', '2.45[GHz]', 'H2 plasma frequency');
model.param.set('n_plasma', '1e17[1/m^3]', 'Plasma density');
model.param.set('T_plasma', '5[eV]', 'Plasma temperature');
model.param.set('V_aln', '100[V]', 'AlN actuation voltage');
model.param.set('f_aln', '100[MHz]', 'AlN actuation frequency');
model.param.set('I_coil', '10[A]', 'Coil current');
model.param.set('lambda_laser', '1.55[um]', 'Femtolaser wavelength');

% Harmonic nodes (odd angles: 30, 90, 150, 210, 270 degrees)
model.param.set('node_1', '30[deg]');
model.param.set('node_3', '90[deg]');
model.param.set('node_5', '150[deg]');
model.param.set('node_7', '210[deg]');
model.param.set('node_9', '270[deg]');

%% Component
model.component.create('comp1', true);

%% Geometry
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').geom('geom1').lengthUnit('um');

% 1. Mo nozzle (cylindrical, z-axis oriented)
model.component('comp1').geom('geom1').create('nozzle', 'Cylinder');
model.component('comp1').geom('geom1').feature('nozzle').set('r', 'r_nozzle_outer');
model.component('comp1').geom('geom1').feature('nozzle').set('h', 'h_nozzle');
model.component('comp1').geom('geom1').feature('nozzle').set('pos', [0 0 0]);

model.component('comp1').geom('geom1').create('nozzle_inner', 'Cylinder');
model.component('comp1').geom('geom1').feature('nozzle_inner').set('r', 'r_nozzle_inner');
model.component('comp1').geom('geom1').feature('nozzle_inner').set('h', 'h_nozzle+10[um]');
model.component('comp1').geom('geom1').feature('nozzle_inner').set('pos', [0 0 -5]);

model.component('comp1').geom('geom1').create('nozzle_diff', 'Difference');
model.component('comp1').geom('geom1').feature('nozzle_diff').selection('input').set('nozzle');
model.component('comp1').geom('geom1').feature('nozzle_diff').selection('input2').set('nozzle_inner');

% 2. CTE Toroid (SiC backbone, CVD middle, c-BN outer)
model.component('comp1').geom('geom1').create('toroid', 'Torus');
model.component('comp1').geom('geom1').feature('toroid').set('rmaj', 'r_toroid_major');
model.component('comp1').geom('geom1').feature('toroid').set('rmin', 'r_toroid_minor');
model.component('comp1').geom('geom1').feature('toroid').set('pos', [0 0 'h_nozzle+gap_height']);

% 3. Spiral Mo coils (helical structure)
model.component('comp1').geom('geom1').create('coil_path', 'Helix');
model.component('comp1').geom('geom1').feature('coil_path').set('r', 'r_coil');
model.component('comp1').geom('geom1').feature('coil_path').set('h', 'h_coil');
model.component('comp1').geom('geom1').feature('coil_path').set('pitch', 'pitch_coil');
model.component('comp1').geom('geom1').feature('coil_path').set('nturns', 25);
model.component('comp1').geom('geom1').feature('coil_path').set('pos', [0 0 'h_nozzle+gap_height+r_toroid_minor']);

% 4. Hollow sapphire resonator
model.component('comp1').geom('geom1').create('resonator_outer', 'Cylinder');
model.component('comp1').geom('geom1').feature('resonator_outer').set('r', 'r_resonator_outer');
model.component('comp1').geom('geom1').feature('resonator_outer').set('h', 'h_resonator');
model.component('comp1').geom('geom1').feature('resonator_outer').set('pos', [0 0 'h_nozzle+gap_height+r_toroid_minor+h_coil']);

model.component('comp1').geom('geom1').create('resonator_inner', 'Cylinder');
model.component('comp1').geom('geom1').feature('resonator_inner').set('r', 'r_resonator_inner');
model.component('comp1').geom('geom1').feature('resonator_inner').set('h', 'h_resonator+10[um]');
model.component('comp1').geom('geom1').feature('resonator_inner').set('pos', [0 0 'h_nozzle+gap_height+r_toroid_minor+h_coil-5[um]']);

model.component('comp1').geom('geom1').create('resonator_diff', 'Difference');
model.component('comp1').geom('geom1').feature('resonator_diff').selection('input').set('resonator_outer');
model.component('comp1').geom('geom1').feature('resonator_diff').selection('input2').set('resonator_inner');

model.component('comp1').geom('geom1').run;

%% Materials
% Molybdenum (Mo) - for nozzle and coils
model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').label('Molybdenum (Mo)');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', '1.89e7[S/m]');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', '138[W/(m*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('density', '10280[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', '251[J/(kg*K)]');

% Aluminum Nitride (AlN) - piezoelectric actuator
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').label('Aluminum Nitride (AlN)');
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', '9.14');
model.component('comp1').material('mat2').propertyGroup('def').set('thermalconductivity', '285[W/(m*K)]');
model.component('comp1').material('mat2').propertyGroup('def').set('density', '3260[kg/m^3]');
model.component('comp1').material('mat2').propertyGroup('def').set('d33', '5.53e-12[m/V]');

% Silicon Carbide (SiC) - toroid backbone
model.component('comp1').material.create('mat3', 'Common');
model.component('comp1').material('mat3').label('Silicon Carbide (SiC)');
model.component('comp1').material('mat3').propertyGroup('def').set('thermalconductivity', '490[W/(m*K)]');
model.component('comp1').material('mat3').propertyGroup('def').set('density', '3210[kg/m^3]');
model.component('comp1').material('mat3').propertyGroup('def').set('relpermittivity', '9.7');

% Cubic Boron Nitride (c-BN) - toroid outer layer
model.component('comp1').material.create('mat4', 'Common');
model.component('comp1').material('mat4').label('Cubic Boron Nitride (c-BN)');
model.component('comp1').material('mat4').propertyGroup('def').set('thermalconductivity', '1300[W/(m*K)]');
model.component('comp1').material('mat4').propertyGroup('def').set('density', '3480[kg/m^3]');

% Sapphire (Al2O3) - resonator
model.component('comp1').material.create('mat5', 'Common');
model.component('comp1').material('mat5').label('Sapphire (Al2O3)');
model.component('comp1').material('mat5').propertyGroup('def').set('relpermittivity', '11.5');
model.component('comp1').material('mat5').propertyGroup('def').set('thermalconductivity', '46[W/(m*K)]');
model.component('comp1').material('mat5').propertyGroup('def').set('density', '3980[kg/m^3]');
model.component('comp1').material('mat5').propertyGroup('def').set('refractiveindex', '1.76');

% H2 Plasma
model.component('comp1').material.create('mat6', 'Common');
model.component('comp1').material('mat6').label('H2 Plasma');
model.component('comp1').material('mat6').propertyGroup('def').set('relpermittivity', '1-5.64e4*n_plasma/(f_plasma^2)');

%% Physics Modules

% 1. Electromagnetic Waves, Frequency Domain (emw)
model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw').prop('ShapeProperty').set('order_electricfield', 2);
model.component('comp1').physics('emw').feature('wee1').set('f', 'f_plasma');

% Port for plasma input at nozzle
model.component('comp1').physics('emw').create('port1', 'Port', 2);
model.component('comp1').physics('emw').feature('port1').set('PortType', 'Rectangular');
model.component('comp1').physics('emw').feature('port1').set('Pport', '1[kW]');

% Scattering boundary condition
model.component('comp1').physics('emw').create('sbc1', 'ScatteringBoundaryCondition', 2);

% 2. Plasma Module (plas)
model.component('comp1').physics.create('plas', 'Plasma', 'geom1');
model.component('comp1').physics('plas').create('pfd1', 'PlasmaFluidDynamics', 3);
model.component('comp1').physics('plas').feature('pfd1').set('ne', 'n_plasma');
model.component('comp1').physics('plas').feature('pfd1').set('Te', 'T_plasma');

% 3. Heat Transfer in Solids (ht)
model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');

% Cryogenic cooling at resonator
model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 2);
model.component('comp1').physics('ht').feature('temp1').set('T0', 'T_cryo');

% Heat source from plasma
model.component('comp1').physics('ht').create('hs1', 'HeatSource', 3);
model.component('comp1').physics('ht').feature('hs1').set('Q0', '1e6[W/m^3]');

% 4. Piezoelectric Effect (pzd)
model.component('comp1').physics.create('pzd', 'Piezoelectricity', 'geom1');
model.component('comp1').physics('pzd').create('term1', 'Terminal', 2);
model.component('comp1').physics('pzd').feature('term1').set('V0', 'V_aln*sin(2*pi*f_aln*t)');

% 5. Magnetic Fields (mf) for coils
model.component('comp1').physics.create('mf', 'MagneticFields', 'geom1');
model.component('comp1').physics('mf').create('coil1', 'Coil', 3);
model.component('comp1').physics('mf').feature('coil1').set('ICoil', 'I_coil*cos(2*pi*f_aln*t)');
model.component('comp1').physics('mf').feature('coil1').set('CoilType', 'Helical');

%% Multiphysics Couplings
model.component('comp1').multiphysics.create('emh1', 'ElectromagneticHeating', 3);
model.component('comp1').multiphysics.create('tee1', 'ThermalExpansion', 3);

%% Mesh
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').autoMeshSize(4); % Fine mesh

%% Studies

% Study 1: Eigenfrequency for spiral plasma modes
model.study.create('std1');
model.study('std1').create('eig', 'Eigenfrequency');
model.study('std1').label('Eigenmode Analysis - H2 Plasma Spiral');
model.study('std1').feature('eig').set('neigs', 20);
model.study('std1').feature('eig').set('shift', 'f_plasma');

% Study 2: Frequency Domain for EM wave propagation
model.study.create('std2');
model.study('std2').create('freq', 'Frequency');
model.study('std2').label('EM Wave Propagation');
model.study('std2').feature('freq').set('plist', 'range(1e9,0.1e9,5e9)');

% Study 3: Time Dependent for polariton dynamics
model.study.create('std3');
model.study('std3').create('time', 'Transient');
model.study('std3').label('Polariton Formation and Avoided Crossings');
model.study('std3').feature('time').set('tlist', 'range(0,1e-12,100e-12)');

% Study 4: Harmonic Analysis on odd nodes (1,3,5,7...270 degrees)
model.study.create('std4');
model.study('std4').create('freq', 'Frequency');
model.study('std4').label('Harmonic Analysis - Odd Nodes');
model.study('std4').feature('freq').set('plist', 'f_aln*[1,3,5,7,9]');

% Study 5: Stationary for Casimir effect in whispering gallery modes
model.study.create('std5');
model.study('std5').create('stat', 'Stationary');
model.study('std5').label('Whispering Gallery Modes - Casimir Effect');

%% Results Visualization

% Poynting vector visualization
model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').label('Poynting Vector - Directional Bias');
model.result('pg1').create('slc1', 'Slice');
model.result('pg1').feature('slc1').set('expr', 'emw.normPoav');
model.result('pg1').create('arw1', 'Arrow3D');
model.result('pg1').feature('arw1').set('expr', {'emw.Poavx', 'emw.Poavy', 'emw.Poavz'});

% Electric field for polariton analysis
model.result.create('pg2', 'PlotGroup3D');
model.result('pg2').label('Polariton Field Distribution');
model.result('pg2').create('slc1', 'Slice');
model.result('pg2').feature('slc1').set('expr', 'emw.normE');

% Temperature distribution
model.result.create('pg3', 'PlotGroup3D');
model.result('pg3').label('Temperature Distribution');
model.result('pg3').create('slc1', 'Slice');
model.result('pg3').feature('slc1').set('expr', 'T');

% Plasma density
model.result.create('pg4', 'PlotGroup3D');
model.result('pg4').label('H2 Plasma Density');
model.result('pg4').create('slc1', 'Slice');
model.result('pg4').feature('slc1').set('expr', 'plas.ne');

% Eigenmode visualization
model.result.create('pg5', 'PlotGroup3D');
model.result('pg5').label('Spiral Eigenmode Shapes');
model.result('pg5').create('slc1', 'Slice');
model.result('pg5').feature('slc1').set('expr', 'abs(emw.Ex)');

% Thrust calculation via Poynting vector
model.result.numerical.create('int1', 'IntSurface');
model.result.numerical('int1').label('Net Thrust from Poynting Vector');
model.result.numerical('int1').set('expr', 'emw.Poavz');
model.result.numerical('int1').set('descr', 'Z-component thrust');

%% Export Configuration
model.result.export.create('data1', 'Data');
model.result.export('data1').label('Export All Results');
model.result.export('data1').set('filename', 'results/comsol_drakarius_data.txt');

%% Display Information
fprintf('\n=== DRAKARIUS COMSOL MODEL CREATED ===\n');
fprintf('Model: %s\n', model.label);
fprintf('Studies configured:\n');
fprintf('  1. Eigenmode Analysis - H2 Plasma Spiral\n');
fprintf('  2. EM Wave Propagation\n');
fprintf('  3. Polariton Formation and Avoided Crossings\n');
fprintf('  4. Harmonic Analysis - Odd Nodes\n');
fprintf('  5. Whispering Gallery Modes - Casimir Effect\n');
fprintf('\nTo run simulation:\n');
fprintf('  model.study(''std1'').run;  %% Eigenmode\n');
fprintf('  model.study(''std2'').run;  %% EM Waves\n');
fprintf('  model.study(''std3'').run;  %% Polaritons\n');
fprintf('  model.study(''std4'').run;  %% Harmonics\n');
fprintf('  model.study(''std5'').run;  %% Casimir/WGM\n');
fprintf('\nTo save:\n');
fprintf('  mphsave(model, ''results/drakarius_full_system.mph'');\n');

end
