# drakarius-validation
Open-source simulations validating the Drakarius Framework for resonant plasma-lattice coupling and crystalline vortex propulsion. See manuscript on Zenodo DOI: 10.5281/zenodo.17289779.

## Overview

This repository contains comprehensive multiphysics simulations for the Drakarius propulsion and energy system, including:

- **Python Simulations**: Level 1 (1D polariton vortex) and Level 3 (3D full propulsion)
- **COMSOL Multiphysics**: Complete system simulation with all physics modules

## COMSOL Multiphysics Simulation

The `simulations/` directory contains a comprehensive COMSOL Multiphysics model implementing:

### System Components
1. AlN-actuated Mo nozzle (z-axis cylindrical) - eigenmoded H2 plasma spirals
2. CTE Toroid (SiC backbone, CVD middle, c-BN outer plasma-facing)
3. AlN-actuated Mo spiral coils (intex/outex modes at 30°, 90°, 150°, 210°, 270°)
4. Hollow sapphire CVD resonator with femtolaser-cut recursive jewel beetle lattices (cooled to 4K)

### Physics Modules
1. **Eigenmode Analysis** - H2 plasma spiral modes in nozzle
2. **Electromagnetic Wave Propagation** - Frequency domain (1-5 GHz)
3. **Plasma Fluid Dynamics** - H2 species, drift-diffusion transport
4. **Heat Transfer** - Cryogenic cooling (4K), plasma heating, EM coupling
5. **Piezoelectric Actuation** - AlN actuators @ 100V, 100 MHz
6. **Magnetic Fields** - Helical coil dynamics (intex/outex)
7. **Harmonic Analysis** - Odd nodes (1,3,5,7,9 harmonics)
8. **Polariton Formation** - Avoided crossings, Rabi splitting
9. **Whispering Gallery Modes** - Casimir effect in resonator
10. **Poynting Vector Analysis** - Directional bias, thrust calculation

### Files
- `comsol_drakarius_full.java` - Java API implementation
- `comsol_drakarius_full.m` - MATLAB LiveLink implementation
- `comsol_setup_and_analysis.py` - Python analysis and visualization
- `COMSOL_README.md` - Detailed documentation

### Quick Start
```bash
# Generate visualizations and configuration
cd simulations
python comsol_setup_and_analysis.py

# Run COMSOL (requires COMSOL license)
comsol batch -inputfile comsol_drakarius_full.java
# OR in MATLAB with LiveLink:
# model = comsol_drakarius_full();
```

## Python Simulations

- `level1_vortex_propulsion.py` - 1D polariton vortex propagation
- `level3_full_propulsion.py` - 3D plasma-polariton-piezo propulsion

## Results

Results are saved in `simulations/results/`:
- COMSOL model: `drakarius_full_system.mph`
- Configuration: `comsol_config.json`
- Summary: `simulation_summary.txt`
- Visualizations: PNG files (geometry, dispersion, harmonics, WGM)

## References

- Drakarius Framework Manuscript: [DOI 10.5281/zenodo.17289779](https://doi.org/10.5281/zenodo.17289779)
- COMSOL Multiphysics Documentation
- Polariton Physics & Avoided Crossings
- Casimir Effect in Microcavities
- Whispering Gallery Modes
