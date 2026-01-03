# COMSOL Multiphysics Simulation for Drakarius Propulsion System

This directory contains comprehensive COMSOL Multiphysics simulation files for the Drakarius propulsion and energy system.

## Overview

The Drakarius system is a novel plasma-polariton-piezoelectric propulsion device consisting of:

1. **AlN-actuated Mo nozzle** (z-axis oriented cylindrical) - Creates eigenmoded H2 plasma spirals
2. **CTE Toroid** (SiC backbone, CVD middle, c-BN outer) - Vortexes polaritons upward
3. **Spiral Mo coils** - Actuated to spiral polaritons up the z-axis (intex/outex modes)
4. **Hollow sapphire resonator** (CVD) - Femtolaser-cut with recursive jewel beetle lattices, cooled to 4K

## Files

### COMSOL Model Files
- `comsol_drakarius_full.java` - Java API implementation (primary)
- `comsol_drakarius_full.m` - MATLAB LiveLink implementation
- `comsol_setup_and_analysis.py` - Python analysis and documentation tool

### Output Files (generated)
- `results/drakarius_full_system.mph` - Complete COMSOL model
- `results/comsol_drakarius_data.txt` - Exported simulation data
- `results/comsol_config.json` - Configuration parameters
- `results/simulation_summary.txt` - Comprehensive summary
- `results/*.png` - Visualization plots

## Physics Modules Implemented

### 1. Eigenmode Analysis
- **Purpose**: Compute H2 plasma spiral eigenmodes in the narrow Mo nozzle
- **Module**: Electromagnetic Waves - Eigenfrequency
- **Key Features**:
  - 20 eigenmodes around plasma frequency (2.45 GHz)
  - Spiral mode formation for z-axis propagation
  - m-fold azimuthal symmetry

### 2. Electromagnetic Wave Propagation
- **Purpose**: Model EM wave dynamics through the system
- **Module**: Electromagnetic Waves - Frequency Domain
- **Key Features**:
  - Frequency sweep: 1-5 GHz
  - Port excitation at nozzle (1 kW input power)
  - Scattering boundary conditions
  - Poynting vector calculation for thrust

### 3. Plasma Fluid Dynamics
- **Purpose**: Simulate H2 plasma behavior
- **Module**: Plasma Physics
- **Key Features**:
  - Species: H2
  - Density: 10¹⁷ m⁻³
  - Temperature: 5 eV
  - Drift-diffusion transport
  - Coupling to EM fields

### 4. Heat Transfer
- **Purpose**: Thermal management and distribution
- **Module**: Heat Transfer in Solids
- **Key Features**:
  - Cryogenic cooling (4K) at sapphire resonator
  - Plasma heating sources
  - Material thermal properties
  - Multiphysics coupling to EM heating

### 5. Piezoelectric Actuation
- **Purpose**: AlN actuator dynamics
- **Module**: Piezoelectricity
- **Key Features**:
  - Sinusoidal voltage: 100V @ 100 MHz
  - AlN material (d₃₃ = 5.53 pC/N)
  - Applied to nozzle, toroid, coils, and resonator
  - Coupling to structural mechanics

### 6. Magnetic Fields
- **Purpose**: Spiral coil electromagnetic dynamics
- **Module**: Magnetic Fields
- **Key Features**:
  - Helical coil geometry
  - Time-varying current: 10A cos(2πft)
  - Intex mode: spiral polaritons upward
  - Outex mode: spiral polaritons downward

### 7. Harmonic Analysis
- **Purpose**: Study odd harmonic nodes
- **Module**: Frequency Domain
- **Key Features**:
  - Nodes at: 30°, 90°, 150°, 210°, 270°
  - Harmonics: 1, 3, 5, 7, 9 × base frequency
  - Resonance identification

### 8. Polariton Formation
- **Purpose**: Light-matter coupling and avoided crossings
- **Module**: Time-Dependent Transient
- **Key Features**:
  - Exciton-photon coupling
  - Upper and lower polariton branches
  - Avoided crossing calculation
  - Rabi splitting analysis

### 9. Whispering Gallery Modes & Casimir Effect
- **Purpose**: Resonator quantum effects
- **Module**: Stationary Analysis
- **Key Features**:
  - Sapphire resonator @ 4K
  - Femtolaser-cut recursive lattices
  - WGM resonances
  - Casimir force in confined geometry
  - Enhanced coherence

### 10. Poynting Vector Analysis
- **Purpose**: Directional bias and thrust calculation
- **Module**: Derived Values
- **Key Features**:
  - 3D Poynting vector field: **S** = **E** × **H**
  - Surface integration for net thrust
  - Z-component momentum flux
  - Validation of propulsion effect

## Multiphysics Couplings

The simulation includes critical multiphysics interactions:

1. **Electromagnetic Heating**: EM field energy → Heat generation
2. **Thermal Expansion**: Temperature → Mechanical stress
3. **Piezoelectric Coupling**: Voltage → Strain → EM field modification
4. **Plasma-EM Coupling**: Plasma density → Relative permittivity

## Material Properties

### Molybdenum (Mo)
- Conductivity: 1.89×10⁷ S/m
- Thermal conductivity: 138 W/(m·K)
- Density: 10,280 kg/m³

### Aluminum Nitride (AlN)
- Relative permittivity: 9.14
- Thermal conductivity: 285 W/(m·K)
- Piezoelectric coefficient d₃₃: 5.53 pm/V

### Silicon Carbide (SiC)
- Thermal conductivity: 490 W/(m·K)
- Relative permittivity: 9.7

### Cubic Boron Nitride (c-BN)
- Thermal conductivity: 1300 W/(m·K) (highest thermal conductivity)
- Density: 3,480 kg/m³

### Sapphire (Al₂O₃)
- Relative permittivity: 11.5
- Refractive index: 1.76
- Thermal conductivity: 46 W/(m·K)

### H2 Plasma
- Permittivity: εᵣ = 1 - ωₚ²/ω²
- Density-dependent dispersion

## Running the Simulation

### Method 1: Java API (Recommended)
```bash
# Requires COMSOL with Java API
comsol batch -inputfile comsol_drakarius_full.java
```

### Method 2: MATLAB LiveLink
```matlab
% Start COMSOL server
% In MATLAB:
model = comsol_drakarius_full();

% Run studies
model.study('std1').run;  % Eigenmode
model.study('std2').run;  % EM Waves
model.study('std3').run;  % Polaritons
model.study('std4').run;  % Harmonics
model.study('std5').run;  % Casimir/WGM

% Save
mphsave(model, 'results/drakarius_full_system.mph');
```

### Method 3: Python Analysis (No COMSOL required)
```bash
python comsol_setup_and_analysis.py
```
This generates:
- Geometry visualization
- Polariton dispersion curves
- Harmonic node analysis
- WGM spectrum
- Configuration files
- Summary report

## Studies Configured

1. **Study 1: Eigenfrequency Analysis**
   - Computes 20 eigenmodes
   - Search near plasma frequency
   - Spiral mode shapes

2. **Study 2: Frequency Domain**
   - Sweep: 1-5 GHz (100 MHz steps)
   - EM wave propagation
   - S-parameters

3. **Study 3: Time-Dependent**
   - Duration: 0-100 ps
   - Time step: 1 ps
   - Polariton dynamics

4. **Study 4: Harmonic Analysis**
   - Odd harmonics: f × [1,3,5,7,9]
   - Node locations: 30°, 90°, 150°, 210°, 270°

5. **Study 5: Stationary**
   - Steady-state WGM
   - Casimir effect calculation

## Expected Results

### Validation Metrics
- ✓ Eigenmode spiral structure with m-fold symmetry
- ✓ Polariton avoided crossing gap: ~ℏΩᵣ (Rabi splitting)
- ✓ Positive z-directed Poynting vector (thrust)
- ✓ Energy conservation: ΔE/E < 1%
- ✓ Temperature stability: T ≈ 4K in resonator
- ✓ Whispering gallery Q-factor > 10⁶

### Thrust Calculation
Net thrust = ∫∫ Sᵧ/c dA (Poynting vector flux)

### Polariton Properties
- Avoided crossing at k₀
- Rabi splitting: Ω ≈ 0.1% of ω₀
- Group velocity: ~0.95c

## Additional Validation Tests

The simulation framework supports additional studies:

1. **Parameter Sweep**: Vary plasma density, actuation voltage, current
2. **Sensitivity Analysis**: Material property variations
3. **Optimization**: Geometry optimization for maximum thrust
4. **Transient Startup**: System initialization dynamics
5. **Thermal Cycling**: 4K to 300K temperature effects

## References

- Drakarius Framework: DOI 10.5281/zenodo.17289779
- COMSOL Multiphysics Documentation
- Polariton Physics: Kavokin et al., "Microcavities"
- Casimir Effect: Milton, "The Casimir Effect"
- Whispering Gallery Modes: Vahala, "Optical Microcavities"

## Notes

- Mesh: Fine (element size ~λ/10 at highest frequency)
- Solver: Direct (MUMPS) for eigenmode, Iterative for frequency
- Memory requirements: ~16-32 GB RAM
- Computation time: ~2-12 hours depending on hardware
- Parallel processing: Recommended for transient study

## Contact

For questions about the Drakarius simulation framework, see:
- Repository: https://github.com/AmenRaAsar/drakarius-validation
- Manuscript: DOI 10.5281/zenodo.17289779
