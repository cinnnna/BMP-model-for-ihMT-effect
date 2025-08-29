# BMP-model-for-ihMT-effect

This file provides guidance to utilise the code to reproduce the BMP model for inhomogeneous magnetization transfer (ihMT) effect. This is a MATLAB codebase implementing a Bloch-McConnell model for simulating ihMT effects in MRI. The code models magnetic resonance interactions between free water pools and semisolid pools using multi-component tissue models.

## Key Files and Architecture

### Core Simulation Engine
- `ihMT_integrate.m` - Main integration function that propagates magnetization through Bloch-McConnell equations during RF pulses
- `init_tissue.m` - Tissue parameter initialization for different tissue types ('hc' for hair conditioner, 'egg' for egg white, 'random' for fitting)
- `Z_spectrum_simulation.m` - Simulates Z-spectra by sweeping frequency offsets

### Pulse Generation and Lineshapes
- `gen_pulse.m` - Generates RF pulse shapes (square/gaussian) with modulation for single-band or dual-band approaches
- `SuperLorentzian_lineshape.m` - Computes Super-Lorentzian lineshape for semisolid pools
- `gauss_lineshape.m` - Computes Gaussian lineshape alternative

### Data Fitting Scripts
- `Fitting_egg_sample.m` - Simultaneous fitting of experimental egg white datasets with parameter bounds
- `Fitting_hc_sample.m` - Fitting script for hair conditioner datasets

### Testing and Validation
- `steadystate_check.m` - Tests steady-state behavior across different RF powers

## Model Architecture

The codebase implements a 4-pool Bloch-McConnell model:
1. **Free water pool** - characterized by R1, R2 rates
2. **Semisolid Zeeman pool** - with M0, R1, T2 parameters
3. **Semisolid dipolar pool** - with R1D (dipolar relaxation) parameter
4. **Exchange processes** - governed by rate constant k

### Key Parameters
- Pulse parameters: duration (typically 5s), B1 amplitude (μT), shape (square/gaussian)
- Tissue parameters: R1f, R2f, M0s, R1s, R1D, f (dipolar fraction), T2s, k
- Frequency offsets: typically swept from 0-6 kHz for Z-spectra

### Band Configurations
- `'1band'` - Single-band RF irradiation with complex modulation
- `'2band'` - Dual-band RF with cosine modulation (√2 scaling factor applied)

## Running the Code

Since this is MATLAB code, execute scripts directly in MATLAB:

```matlab
% Run Z-spectrum simulation
Z_spectrum_simulation

% Test steady-state behavior
steadystate_check

% Fit experimental data
Fitting_egg_sample
```

No build system, package management, or external dependencies beyond standard MATLAB functions are used.

## Development Notes

- All magnetic field units are in \(\mu\)T (microtesla)
- Time units are in seconds
- Frequency offsets are in Hz
- Relaxation rates are in \(s^{-1}\)
- The model uses matrix exponentiation (`expm`) for numerical integration of differential equations
- Parameter bounds and initial values are defined within fitting scripts
