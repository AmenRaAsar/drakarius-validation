"""
COMSOL Multiphysics Model Documentation for Drakarius Propulsion System
This Python script documents the complete physics setup and provides analysis tools

Physics Modules Implemented:
1. Eigenmode analysis for H2 plasma spiral modes
2. Electromagnetic wave propagation
3. Plasma fluid dynamics
4. Heat transfer with cryogenic cooling (4K)
5. Piezoelectric actuation (AlN)
6. Harmonic studies on odd nodes (1,3,5,7...270°)
7. Casimir effect in whispering gallery modes
8. Polariton formation and avoided crossings
9. Poynting vector directional bias analysis

Author: Drakarius Framework Validation Team
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from mpl_toolkits.mplot3d import Axes3D
import json

# Physical constants
c0 = 299792458  # Speed of light (m/s)
eps0 = 8.854187817e-12  # Permittivity of free space (F/m)
mu0 = 1.256637061e-6  # Permeability of free space (H/m)
e_charge = 1.60217662e-19  # Elementary charge (C)
h_planck = 6.62607015e-34  # Planck constant (J·s)
hbar = h_planck / (2 * np.pi)  # Reduced Planck constant
k_B = 1.380649e-23  # Boltzmann constant (J/K)

class DrakariusSimulationConfig:
    """Configuration class for Drakarius COMSOL simulation"""
    
    def __init__(self):
        # Geometric parameters (meters)
        self.r_nozzle_inner = 50e-6  # Mo nozzle inner radius
        self.r_nozzle_outer = 100e-6  # Mo nozzle outer radius
        self.h_nozzle = 500e-6  # Mo nozzle height
        self.gap_height = 200e-6  # Gap to toroid
        self.r_toroid_major = 2e-3  # CTE toroid major radius
        self.r_toroid_minor = 300e-6  # CTE toroid minor radius
        self.h_coil = 5e-3  # Height of spiral coils
        self.r_coil = 1.5e-3  # Radius of spiral coils
        self.pitch_coil = 200e-6  # Coil pitch
        self.r_resonator_outer = 3e-3  # Sapphire resonator outer radius
        self.r_resonator_inner = 2.5e-3  # Sapphire resonator inner radius
        self.h_resonator = 4e-3  # Resonator height
        
        # Physics parameters
        self.T_cryo = 4  # Cryogenic temperature (K)
        self.f_plasma = 2.45e9  # H2 plasma frequency (Hz)
        self.n_plasma = 1e17  # Plasma density (1/m^3)
        self.T_plasma_eV = 5  # Plasma temperature (eV)
        self.T_plasma_K = self.T_plasma_eV * 11600  # Plasma temperature (K)
        self.V_aln = 100  # AlN actuation voltage (V)
        self.f_aln = 100e6  # AlN actuation frequency (Hz)
        self.I_coil = 10  # Coil current (A)
        self.lambda_laser = 1.55e-6  # Femtolaser wavelength (m)
        
        # Material properties
        self.materials = {
            'Mo': {
                'conductivity': 1.89e7,  # S/m
                'thermal_conductivity': 138,  # W/(m·K)
                'density': 10280,  # kg/m^3
                'heat_capacity': 251  # J/(kg·K)
            },
            'AlN': {
                'permittivity': 9.14,
                'thermal_conductivity': 285,  # W/(m·K)
                'density': 3260,  # kg/m^3
                'd33': 5.53e-12  # Piezoelectric coefficient (m/V)
            },
            'SiC': {
                'thermal_conductivity': 490,  # W/(m·K)
                'density': 3210,  # kg/m^3
                'permittivity': 9.7
            },
            'c-BN': {
                'thermal_conductivity': 1300,  # W/(m·K)
                'density': 3480  # kg/m^3
            },
            'Sapphire': {
                'permittivity': 11.5,
                'thermal_conductivity': 46,  # W/(m·K)
                'density': 3980,  # kg/m^3
                'refractive_index': 1.76
            }
        }
        
        # Harmonic nodes (odd angles in degrees)
        self.harmonic_nodes = [30, 90, 150, 210, 270]
        
    def calculate_plasma_frequency(self):
        """Calculate plasma frequency from density"""
        omega_p = np.sqrt(self.n_plasma * e_charge**2 / (eps0 * 9.1093837e-31))
        return omega_p / (2 * np.pi)
    
    def calculate_debye_length(self):
        """Calculate Debye length"""
        lambda_D = np.sqrt(eps0 * k_B * self.T_plasma_K / (self.n_plasma * e_charge**2))
        return lambda_D
    
    def calculate_polariton_dispersion(self, k_values):
        """
        Calculate polariton dispersion relation for avoided crossings
        k_values: array of wave vectors (1/m)
        Returns: upper and lower polariton branch frequencies
        """
        omega_c = 2 * np.pi * c0 / self.lambda_laser
        omega_ex = omega_c  # Exciton frequency
        g_coupling = 0.001 * omega_c  # Rabi splitting
        
        # Photon dispersion
        omega_ph = c0 * k_values / np.sqrt(self.materials['Sapphire']['permittivity'])
        
        # Coupled oscillator model
        omega_plus = 0.5 * (omega_ph + omega_ex + np.sqrt((omega_ph - omega_ex)**2 + 4*g_coupling**2))
        omega_minus = 0.5 * (omega_ph + omega_ex - np.sqrt((omega_ph - omega_ex)**2 + 4*g_coupling**2))
        
        return omega_plus / (2*np.pi), omega_minus / (2*np.pi)
    
    def calculate_casimir_force(self):
        """
        Estimate Casimir force between resonator walls
        Using parallel plate approximation
        """
        d = self.r_resonator_outer - self.r_resonator_inner  # Separation
        A = 2 * np.pi * self.r_resonator_inner * self.h_resonator  # Area
        F_casimir = -np.pi**2 * hbar * c0 * A / (240 * d**4)
        return F_casimir
    
    def calculate_whispering_gallery_modes(self, m_max=20):
        """
        Calculate whispering gallery mode frequencies
        m_max: maximum azimuthal mode number
        Returns: array of WGM frequencies
        """
        n_eff = self.materials['Sapphire']['refractive_index']
        r = (self.r_resonator_outer + self.r_resonator_inner) / 2
        
        frequencies = []
        for m in range(1, m_max + 1):
            # Approximate WGM frequency
            f_wgm = m * c0 / (2 * np.pi * r * n_eff)
            frequencies.append(f_wgm)
        
        return np.array(frequencies)
    
    def save_config(self, filename='comsol_config.json'):
        """Save configuration to JSON file"""
        config_dict = {
            'geometry': {
                'r_nozzle_inner': self.r_nozzle_inner,
                'r_nozzle_outer': self.r_nozzle_outer,
                'h_nozzle': self.h_nozzle,
                'gap_height': self.gap_height,
                'r_toroid_major': self.r_toroid_major,
                'r_toroid_minor': self.r_toroid_minor,
                'h_coil': self.h_coil,
                'r_coil': self.r_coil,
                'pitch_coil': self.pitch_coil,
                'r_resonator_outer': self.r_resonator_outer,
                'r_resonator_inner': self.r_resonator_inner,
                'h_resonator': self.h_resonator
            },
            'physics': {
                'T_cryo': self.T_cryo,
                'f_plasma': self.f_plasma,
                'n_plasma': self.n_plasma,
                'T_plasma_eV': self.T_plasma_eV,
                'V_aln': self.V_aln,
                'f_aln': self.f_aln,
                'I_coil': self.I_coil,
                'lambda_laser': self.lambda_laser
            },
            'materials': self.materials,
            'harmonic_nodes': self.harmonic_nodes
        }
        
        with open(filename, 'w') as f:
            json.dump(config_dict, f, indent=2)
        
        print(f"Configuration saved to {filename}")


def visualize_geometry(config):
    """Create 2D schematic of Drakarius geometry"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 12))
    
    # Mo nozzle
    nozzle_outer = Rectangle((-config.r_nozzle_outer*1e6, 0), 
                              2*config.r_nozzle_outer*1e6, 
                              config.h_nozzle*1e6,
                              fill=False, edgecolor='blue', linewidth=2, label='Mo Nozzle')
    ax.add_patch(nozzle_outer)
    
    # Gap
    gap_y = config.h_nozzle*1e6
    ax.plot([-config.r_toroid_major*1e3, config.r_toroid_major*1e3], 
            [gap_y, gap_y], 'k--', alpha=0.3, label='Gap')
    
    # Toroid
    toroid_y = gap_y + config.gap_height*1e6
    toroid = Circle((0, toroid_y), config.r_toroid_minor*1e6, 
                    fill=False, edgecolor='green', linewidth=2, label='CTE Toroid')
    ax.add_patch(toroid)
    
    # Coils (simplified as vertical lines)
    coil_y_start = toroid_y + config.r_toroid_minor*1e6
    ax.plot([-config.r_coil*1e3, -config.r_coil*1e3], 
            [coil_y_start, coil_y_start + config.h_coil*1e3], 
            'r-', linewidth=3, label='Spiral Coils')
    ax.plot([config.r_coil*1e3, config.r_coil*1e3], 
            [coil_y_start, coil_y_start + config.h_coil*1e3], 
            'r-', linewidth=3)
    
    # Resonator
    res_y = coil_y_start + config.h_coil*1e3
    res_outer = Rectangle((-config.r_resonator_outer*1e3, res_y), 
                          2*config.r_resonator_outer*1e3, 
                          config.h_resonator*1e3,
                          fill=False, edgecolor='purple', linewidth=2, label='Sapphire Resonator')
    ax.add_patch(res_outer)
    
    res_inner = Rectangle((-config.r_resonator_inner*1e3, res_y), 
                          2*config.r_resonator_inner*1e3, 
                          config.h_resonator*1e3,
                          fill=False, edgecolor='purple', linewidth=1, linestyle='--')
    ax.add_patch(res_inner)
    
    ax.set_xlim(-4, 4)
    ax.set_ylim(-0.5, res_y + config.h_resonator*1e3 + 1)
    ax.set_xlabel('Radial Distance (mm)', fontsize=12)
    ax.set_ylabel('Z-axis (µm to mm)', fontsize=12)
    ax.set_title('Drakarius Propulsion System - Geometry', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig('results/drakarius_geometry.png', dpi=300)
    print("Geometry visualization saved to results/drakarius_geometry.png")
    plt.close()


def visualize_polariton_dispersion(config):
    """Visualize polariton dispersion with avoided crossing"""
    k_values = np.linspace(1e6, 1e8, 1000)
    omega_upper, omega_lower = config.calculate_polariton_dispersion(k_values)
    
    plt.figure(figsize=(10, 6))
    plt.plot(k_values/1e6, omega_upper/1e12, 'b-', linewidth=2, label='Upper Polariton Branch')
    plt.plot(k_values/1e6, omega_lower/1e12, 'r-', linewidth=2, label='Lower Polariton Branch')
    
    # Mark avoided crossing region
    omega_c = c0 / config.lambda_laser
    plt.axhline(omega_c/1e12, color='k', linestyle='--', alpha=0.5, label='Cavity Frequency')
    
    plt.xlabel('Wave Vector k (µm⁻¹)', fontsize=12)
    plt.ylabel('Frequency (THz)', fontsize=12)
    plt.title('Polariton Dispersion Relation - Avoided Crossing', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('results/polariton_dispersion.png', dpi=300)
    print("Polariton dispersion saved to results/polariton_dispersion.png")
    plt.close()


def visualize_harmonic_nodes(config):
    """Visualize harmonic analysis nodes"""
    angles = np.array(config.harmonic_nodes)
    harmonics = np.array([1, 3, 5, 7, 9])
    
    fig = plt.figure(figsize=(14, 6))
    
    # Angular distribution (polar plot)
    ax1 = fig.add_subplot(121, projection='polar')
    theta = np.deg2rad(angles)
    r = harmonics
    ax1.plot(theta, r, 'ro-', markersize=10, linewidth=2)
    ax1.set_theta_zero_location('N')
    ax1.set_title('Harmonic Nodes (Odd)', fontsize=12, fontweight='bold', pad=20)
    
    # Regular axis for frequency spectrum
    ax2 = fig.add_subplot(122)
    
    # Frequency spectrum
    freq_harmonics = config.f_aln * harmonics
    ax2.stem(harmonics, freq_harmonics/1e6, basefmt=' ')
    ax2.set_xlabel('Harmonic Number', fontsize=12)
    ax2.set_ylabel('Frequency (MHz)', fontsize=12)
    ax2.set_title('Harmonic Frequency Spectrum', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('results/harmonic_analysis.png', dpi=300)
    print("Harmonic analysis saved to results/harmonic_analysis.png")
    plt.close()


def visualize_whispering_gallery_modes(config):
    """Visualize whispering gallery mode spectrum"""
    wgm_freqs = config.calculate_whispering_gallery_modes(m_max=30)
    m_values = np.arange(1, len(wgm_freqs) + 1)
    
    plt.figure(figsize=(12, 6))
    plt.plot(m_values, wgm_freqs/1e12, 'o-', markersize=6, linewidth=1.5)
    plt.xlabel('Azimuthal Mode Number m', fontsize=12)
    plt.ylabel('Frequency (THz)', fontsize=12)
    plt.title('Whispering Gallery Modes in Sapphire Resonator (4K)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('results/whispering_gallery_modes.png', dpi=300)
    print("WGM spectrum saved to results/whispering_gallery_modes.png")
    plt.close()


def generate_simulation_summary(config):
    """Generate a comprehensive summary of the simulation"""
    summary = f"""
{'='*80}
DRAKARIUS COMSOL MULTIPHYSICS SIMULATION SUMMARY
{'='*80}

SYSTEM CONFIGURATION:
--------------------
Geometry:
  - Mo Nozzle: Inner radius = {config.r_nozzle_inner*1e6:.1f} µm, Outer = {config.r_nozzle_outer*1e6:.1f} µm
  - Height: {config.h_nozzle*1e6:.1f} µm
  - Gap to Toroid: {config.gap_height*1e6:.1f} µm
  - CTE Toroid: Major radius = {config.r_toroid_major*1e3:.1f} mm, Minor = {config.r_toroid_minor*1e6:.1f} µm
  - Spiral Coils: Height = {config.h_coil*1e3:.1f} mm, Radius = {config.r_coil*1e3:.1f} mm, Pitch = {config.pitch_coil*1e6:.1f} µm
  - Sapphire Resonator: Outer = {config.r_resonator_outer*1e3:.1f} mm, Inner = {config.r_resonator_inner*1e3:.1f} mm
  - Resonator Height: {config.h_resonator*1e3:.1f} mm

Physics Parameters:
  - Cryogenic Temperature: {config.T_cryo} K
  - H2 Plasma Frequency: {config.f_plasma/1e9:.2f} GHz
  - Plasma Density: {config.n_plasma:.2e} m⁻³
  - Plasma Temperature: {config.T_plasma_eV} eV ({config.T_plasma_K:.0f} K)
  - AlN Actuation: {config.V_aln} V @ {config.f_aln/1e6:.0f} MHz
  - Coil Current: {config.I_coil} A
  - Femtolaser Wavelength: {config.lambda_laser*1e6:.2f} µm

Calculated Plasma Properties:
  - Plasma Frequency: {config.calculate_plasma_frequency()/1e9:.2f} GHz
  - Debye Length: {config.calculate_debye_length()*1e6:.2f} µm

Casimir Effect:
  - Gap: {(config.r_resonator_outer - config.r_resonator_inner)*1e6:.0f} µm
  - Estimated Casimir Force: {config.calculate_casimir_force():.2e} N

PHYSICS MODULES IMPLEMENTED:
---------------------------
1. Eigenmode Study
   - H2 plasma spiral eigenmodes in narrow nozzle
   - 20 eigenmodes computed around plasma frequency
   - Spiral mode formation for z-axis propagation

2. Electromagnetic Wave Propagation
   - Frequency domain analysis (1-5 GHz)
   - Port excitation at nozzle base (1 kW)
   - Scattering boundary conditions
   - Poynting vector for thrust calculation

3. Plasma Fluid Dynamics
   - H2 plasma species
   - Drift-diffusion transport
   - Eigenmode coupling

4. Heat Transfer
   - Cryogenic cooling (4K) at resonator
   - Plasma heating source
   - Multimaterial thermal coupling
   - Electromagnetic heating coupling

5. Piezoelectric Actuation (AlN)
   - Sinusoidal voltage excitation
   - Coupling to structural mechanics
   - Actuation at nozzle, toroid, coils, and resonator

6. Magnetic Fields
   - Helical coil configuration
   - Time-varying current (intex/outex modes)
   - Coupling to EM waves

7. Harmonic Analysis
   - Odd harmonic nodes: {', '.join(map(str, config.harmonic_nodes))}°
   - Frequencies: f_aln × [1, 3, 5, 7, 9]

8. Polariton Formation
   - Avoided crossing calculation
   - Upper and lower polariton branches
   - Strong light-matter coupling regime

9. Whispering Gallery Modes
   - Sapphire resonator @ 4K
   - Femtolaser-cut recursive jewel beetle lattices
   - Casimir effect in confined geometry
   - Mode coherence optimization

10. Poynting Vector Analysis
    - Directional bias calculation
    - Net thrust from EM momentum flux
    - Integration over exit plane

MULTIPHYSICS COUPLINGS:
----------------------
- Electromagnetic Heating (EM → Heat)
- Thermal Expansion (Heat → Mechanics)
- Piezoelectric Coupling (Voltage → Mechanics → EM)
- Plasma-EM Coupling (Plasma density → Permittivity)

STUDIES CONFIGURED:
-----------------
Study 1: Eigenfrequency Analysis
Study 2: Frequency Domain (1-5 GHz)
Study 3: Time-Dependent Transient (0-100 ps)
Study 4: Harmonic Analysis (Odd nodes)
Study 5: Stationary (Casimir/WGM)

VALIDATION METRICS:
-----------------
- Eigenmode spiral structure (m-fold symmetry)
- Polariton avoided crossing gap
- Poynting vector z-component (thrust)
- Energy conservation (EM + thermal)
- Temperature stability at 4K
- Whispering gallery Q-factor

OUTPUT FILES:
------------
- COMSOL Model: results/drakarius_full_system.mph
- Data Export: results/comsol_drakarius_data.txt
- Visualizations: results/*.png

REFERENCES:
----------
- Drakarius Framework Manuscript: DOI 10.5281/zenodo.17289779
- AlN Piezoelectric Properties
- Sapphire Optical Properties @ Cryogenic Temperatures
- Casimir Effect in Microcavities
- Polariton Physics in Semiconductor Microcavities

{'='*80}
    """
    
    return summary


def main():
    """Main execution function"""
    print("\n" + "="*80)
    print("DRAKARIUS COMSOL MULTIPHYSICS SIMULATION")
    print("Comprehensive Physics Setup and Analysis")
    print("="*80 + "\n")
    
    # Create configuration
    config = DrakariusSimulationConfig()
    
    # Save configuration
    config.save_config('results/comsol_config.json')
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    visualize_geometry(config)
    visualize_polariton_dispersion(config)
    visualize_harmonic_nodes(config)
    visualize_whispering_gallery_modes(config)
    
    # Generate summary
    summary = generate_simulation_summary(config)
    print(summary)
    
    # Save summary to file
    with open('results/simulation_summary.txt', 'w') as f:
        f.write(summary)
    print("Summary saved to results/simulation_summary.txt")
    
    print("\n" + "="*80)
    print("SETUP COMPLETE")
    print("="*80)
    print("\nNext steps:")
    print("1. Open COMSOL Multiphysics")
    print("2. Run Java script: comsol_drakarius_full.java")
    print("   OR")
    print("3. Run MATLAB script: comsol_drakarius_full.m")
    print("\nThe scripts will:")
    print("  - Build complete 3D geometry")
    print("  - Set up all physics modules")
    print("  - Configure 5 comprehensive studies")
    print("  - Compute eigenmode, frequency, and transient solutions")
    print("  - Export results for validation")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
