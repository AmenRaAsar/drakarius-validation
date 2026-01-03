# COMSOL Simulation Usage Notes

## Important Notes for Users

### Running the Simulations

The COMSOL scripts in this repository require:
- COMSOL Multiphysics 6.0 or later
- Required modules: AC/DC, Heat Transfer, Plasma, RF, Structural Mechanics
- Recommended: 16+ GB RAM for fine mesh
- Computation time: 2-12 hours depending on hardware

### Script Compatibility

1. **Java API (comsol_drakarius_full.java)**
   - Most portable, works with COMSOL batch mode
   - Run with: `comsol batch -inputfile comsol_drakarius_full.java`

2. **MATLAB LiveLink (comsol_drakarius_full.m)**
   - Requires COMSOL LiveLink for MATLAB
   - Interactive MATLAB environment
   - Easiest for debugging and result visualization

3. **Python Analysis (comsol_setup_and_analysis.py)**
   - Does NOT require COMSOL
   - Generates visualizations and configuration
   - Run first to understand the system

### Advanced Validation Scripts

The scripts in `results/` directory are templates that require:
- Active COMSOL model loaded
- Completed initial studies (std1-std5)
- Proper dataset names matching your model

**Important**: These are template scripts to demonstrate the validation approach. They may need modification for your specific COMSOL setup:
- Dataset names (dset1, dset2, etc.) depend on study order
- Parameter names must match your model exactly
- Study names must be created first

### Recommended Workflow

1. **Start with Python analysis**
   ```bash
   python comsol_setup_and_analysis.py
   ```
   Review generated visualizations and config.

2. **Run core COMSOL model**
   - Use Java or MATLAB script to build model
   - Run all 5 core studies (std1-std5)
   - Save the .mph file

3. **Perform validation tests**
   - Use advanced validation scripts as templates
   - Modify dataset names as needed
   - Run parameter sweeps individually

4. **Analyze results**
   - Export data from COMSOL
   - Use Python for post-processing
   - Compare with expected physics

### Known Considerations

1. **COMSOL API Variations**: Dataset naming and parameter setting syntax may vary slightly between COMSOL versions. Consult your COMSOL documentation.

2. **Computational Resources**: The full 3D simulation with fine mesh is computationally intensive. Start with coarse mesh for testing.

3. **Physics Module Requirements**: Ensure all required physics modules are licensed. The simulation uses:
   - Electromagnetic Waves, Frequency Domain
   - Plasma Module
   - Heat Transfer
   - Structural Mechanics (Piezoelectricity)
   - AC/DC Module (Magnetic Fields)

4. **Mesh Requirements**: For accurate results:
   - Element size: λ/10 at highest frequency (~6 mm at 5 GHz)
   - Boundary layer mesh at interfaces
   - Fine mesh in resonator for WGM

### Validation Criteria

Expected results:
- Eigenmode spiral structure with m-fold symmetry
- Polariton avoided crossing gap ~0.1% of ω₀
- Positive z-directed Poynting vector (thrust)
- Energy conservation: ΔE/E < 1%
- Resonator temperature stable at 4K
- Whispering gallery Q-factor > 10⁶

### Troubleshooting

**Memory issues**: Reduce mesh density or use iterative solver
**Convergence problems**: Use smaller time steps or frequency steps
**Long computation times**: Use parametric sweep with coarser parameter grid first
**API errors**: Check COMSOL version compatibility and dataset names

### References

- COMSOL Multiphysics User Guide
- COMSOL API Reference (Java/MATLAB)
- Drakarius Framework: DOI 10.5281/zenodo.17289779

### Support

For COMSOL-specific questions, consult COMSOL support.
For Drakarius physics questions, see the manuscript and repository issues.

---

**Note**: The scripts provided are comprehensive templates demonstrating the complete physics setup. They serve as both working code (for the core model) and documentation (for the validation approaches). Users should adapt them to their specific COMSOL environment and computational resources.
