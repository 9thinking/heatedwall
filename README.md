# Droplet Impact Simulation on Heated Wall

This project contains a Basilisk simulation of a liquid droplet impacting on a heated solid surface. The simulation includes two-phase flow dynamics, surface tension effects, and provides both quantitative and qualitative outputs for analysis.

## Overview

The simulation models:
- **Two-phase flow**: Liquid droplet and surrounding gas
- **Surface tension**: Captures interface dynamics
- **Impact dynamics**: Spreading, splashing, and rebound behavior
- **Heat transfer**: Heated wall boundary condition (framework for thermal coupling)

## Physical Parameters

Default parameters (water-like droplet):
- Droplet diameter: D₀ = 2 mm
- Impact velocity: U₀ = 1.0 m/s
- Liquid density: ρₗ = 1000 kg/m³
- Gas density: ρₘ = 1.0 kg/m³
- Liquid viscosity: μₗ = 1×10⁻³ Pa·s
- Surface tension: σ = 0.072 N/m
- Droplet temperature: T = 293.15 K (20°C)
- Wall temperature: Tᵥᵥₐₗₗ = 373.15 K (100°C)

### Dimensionless Numbers

The simulation automatically calculates:
- **Reynolds number**: Re = ρₗU₀D₀/μₗ (ratio of inertial to viscous forces)
- **Weber number**: We = ρₗU₀²D₀/σ (ratio of inertial to surface tension forces)
- **Ohnesorge number**: Oh = √We/Re (ratio of viscous to inertial-capillary forces)
- **Froude number**: Fr = U₀/√(gD₀) (ratio of inertial to gravitational forces)

## Requirements

### Software Dependencies
- **Basilisk C**: The main simulation framework
  - Installation instructions: http://basilisk.fr/src/INSTALL
  - Requires: gcc, make, and standard Unix utilities
- **gnuplot** (optional): For plotting results
- **ffmpeg** (optional): For generating movies from simulation

### System Requirements
- Linux or macOS (or WSL on Windows)
- 2+ GB RAM recommended
- GCC compiler with C99 support

## Installation

1. **Install Basilisk**:
   ```bash
   # Follow instructions at http://basilisk.fr/src/INSTALL
   cd ~
   git clone https://github.com/basilisk-fr/basilisk.git
   cd basilisk/src
   export BASILISK=$PWD
   ```

2. **Clone this repository**:
   ```bash
   git clone https://github.com/9thinking/heatedwall.git
   cd heatedwall
   ```

3. **Compile the simulation**:
   ```bash
   make
   ```

## Usage

### Basic Execution

Run the simulation with default parameters:
```bash
./droplet_impact
```

Or use the Makefile:
```bash
make run
```

### Custom Parameters

The simulation accepts command-line arguments:
```bash
./droplet_impact [LEVEL] [MAXLEVEL] [U0]
```

Parameters:
- `LEVEL`: Base grid refinement level (default: 9, i.e., 512×512 cells)
- `MAXLEVEL`: Maximum adaptive refinement level (default: 11)
- `U0`: Impact velocity in m/s (default: 1.0)

**Examples**:
```bash
# High-resolution simulation
./droplet_impact 10 12

# Faster impact (2 m/s)
./droplet_impact 9 11 2.0

# Lower resolution for quick testing
./droplet_impact 7 9 1.0
```

### Using Makefile Targets

```bash
make run          # Run with default parameters
make run-hr       # Run with higher resolution (level 10, maxlevel 12)
make run-fast     # Run with 2 m/s impact velocity
make clean        # Remove compiled files
make cleanall     # Remove all output files
make help         # Show help message
```

## Output Files

The simulation generates several output files:

### Quantitative Data

**`quantitative_output.dat`**: Time series data including:
- `t`: Non-dimensional time (in units of D₀/U₀)
- `D/D0`: Spreading factor (instantaneous diameter / initial diameter)
- `xc, yc`: Center of mass position
- `volume`: Droplet volume
- `KE`: Kinetic energy

### Qualitative/Visualization Data

**`interface-XXXX.dat`**: Interface shape at different time steps
- Used for plotting droplet contour
- Output every 0.1 time units

**`field-XXXX.dat`**: Full field data (volume fraction, velocity, pressure)
- Complete field information for detailed analysis
- Can be visualized with Basilisk's built-in tools

**`movie.mp4`**: Animation of the simulation
- Shows volume fraction and velocity field
- Generated automatically during simulation

## Analysis and Visualization

### Quick Data Analysis

View the spreading factor evolution:
```bash
gnuplot -e "set terminal dumb; plot 'quantitative_output.dat' u 1:2 w l title 'D/D0'"
```

### Plotting with Gnuplot

Create a plotting script `plot_results.gp`:
```gnuplot
set terminal png size 800,600
set output 'spreading_factor.png'
set xlabel 'Time (D₀/U₀)'
set ylabel 'Spreading Factor (D/D₀)'
set grid
plot 'quantitative_output.dat' u 1:2 w l lw 2 title 'Spreading'
```

Run: `gnuplot plot_results.gp`

### Understanding the Results

**Key metrics to analyze**:

1. **Maximum spreading factor** (D_max/D₀):
   - Indicates how much the droplet spreads upon impact
   - Typically ranges from 1.5 to 4.0 depending on We and Re

2. **Contact time**:
   - Duration of droplet-surface contact
   - Can be estimated from the time when the droplet leaves the surface

3. **Spreading dynamics**:
   - Initial spreading phase (kinetic energy dominated)
   - Maximum spreading (balance point)
   - Retraction phase (surface tension dominated)
   - Possible rebound or equilibrium

## Physical Insights

### Regime Map

The droplet behavior depends on the Weber and Reynolds numbers:

- **Low We (<5)**: Deposition, minimal spreading
- **Medium We (5-30)**: Spreading without splashing
- **High We (>30)**: Possible splashing and fragmentation
- **High Re**: Inertia-dominated, faster spreading
- **Low Re**: Viscosity-dominated, slower dynamics

### Heat Transfer (Framework)

The current simulation includes temperature parameters for the droplet and wall but does not yet include full thermal coupling. To extend this simulation with heat transfer:

1. Include `#include "diffusion.h"` for heat equation
2. Add temperature scalar field
3. Implement conjugate heat transfer at the wall
4. Calculate Nusselt number and heat flux

## Customization

### Modifying Physical Properties

Edit the parameters in `droplet_impact.c`:
```c
double D0 = 2e-3;           // Droplet diameter
double U0 = 1.0;            // Impact velocity
double RHO_L = 1000.0;      // Liquid density
double MU_L = 1e-3;         // Liquid viscosity
double SIGMA = 0.072;       // Surface tension
```

### Changing Domain Size

Modify in the `init` event:
```c
size(6.0);              // Domain size in units of D0
origin(0., -1.0);       // Origin location
```

### Adjusting Output Frequency

Change the output intervals:
```c
double DT_OUTPUT = 0.1;     // Quantitative output interval
event movies(t = 0; t <= t_end; t += 0.02) // Movie frame rate
```

## Troubleshooting

### Compilation Errors

1. **"qcc: command not found"**:
   - Ensure Basilisk is properly installed
   - Check that `$BASILISK` environment variable is set
   - Add Basilisk to PATH: `export PATH=$PATH:$BASILISK/../bin`

2. **Header file not found**:
   - Verify Basilisk installation directory in Makefile
   - Update `BASILISK` variable if needed

### Runtime Issues

1. **Simulation is very slow**:
   - Reduce grid level: `./droplet_impact 7 9`
   - Decrease `t_end` for shorter simulation time
   - Reduce output frequency

2. **Numerical instability**:
   - Ensure We and Re are in reasonable ranges
   - Increase grid resolution for high We numbers
   - Check that fluid property ratios are not extreme

3. **Movie not generated**:
   - Install ffmpeg: `sudo apt-get install ffmpeg`
   - Or disable movie generation by commenting out the `movies` event

## References

### Basilisk Documentation
- Main site: http://basilisk.fr/
- Tutorial: http://basilisk.fr/Tutorial
- Two-phase flows: http://basilisk.fr/src/two-phase.h

### Scientific Background
- Josserand, C., & Thoroddsen, S. T. (2016). Drop impact on a solid surface. *Annual Review of Fluid Mechanics*, 48, 365-391.
- Yarin, A. L. (2006). Drop impact dynamics: splashing, spreading, receding, bouncing…. *Annual Review of Fluid Mechanics*, 38, 159-192.

## License

This code is provided as-is for educational and research purposes.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Contact

For questions or suggestions, please open an issue on GitHub.
