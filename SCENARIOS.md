# Example Simulation Scenarios for Droplet Impact
#
# This file contains parameters for various droplet impact scenarios
# that can be studied using the simulation.

## Scenario 1: Water droplet at low impact velocity (gentle impact)
# Parameters:
# - D0 = 2 mm
# - U0 = 0.5 m/s
# - Re ≈ 1000, We ≈ 7
# Expected: Spreading without splashing, moderate deformation
# Run: Edit droplet_impact.c and set U0 = 0.5, then compile and run

## Scenario 2: Water droplet at moderate velocity (standard case)
# Parameters:
# - D0 = 2 mm
# - U0 = 1.0 m/s (DEFAULT)
# - Re ≈ 2000, We ≈ 28
# Expected: Significant spreading, possible corona formation
# Run: ./droplet_impact (default parameters)

## Scenario 3: Water droplet at high velocity (splashing regime)
# Parameters:
# - D0 = 2 mm
# - U0 = 2.0 m/s
# - Re ≈ 4000, We ≈ 111
# Expected: Splashing, possible secondary droplet formation
# Run: ./droplet_impact 9 11 2.0

## Scenario 4: High-resolution simulation for detailed physics
# Parameters:
# - Base grid level = 10 (1024 cells)
# - Max level = 12 (adaptive to 4096)
# - Standard velocity
# Expected: More accurate interface resolution, longer computation time
# Run: ./droplet_impact 10 12 1.0

## Scenario 5: Quick test simulation (low resolution)
# Parameters:
# - Base grid level = 7 (128 cells)
# - Max level = 9 (adaptive to 512)
# Expected: Fast results, less accurate
# Run: ./droplet_impact 7 9 1.0

## Physical Property Variations (requires code modification)

### Viscous droplet (e.g., glycerol mixture)
# Edit in droplet_impact.c:
# MU_L = 5e-3;  // 5x water viscosity
# Expected: Lower Re, slower spreading, more viscous dissipation

### Low surface tension (e.g., with surfactant)
# Edit in droplet_impact.c:
# SIGMA = 0.035;  // Half of water
# Expected: Higher We, more spreading, easier splashing

### Larger droplet
# Edit in droplet_impact.c:
# D0 = 5e-3;  // 5 mm diameter
# Expected: Higher Re and We, more complex dynamics

### Smaller droplet
# Edit in droplet_impact.c:
# D0 = 1e-3;  // 1 mm diameter
# Expected: Lower Re and We, surface tension dominated

## Heat Transfer Studies (framework for future extension)

### Cold wall (T_wall < T_droplet)
# T_WALL = 273.15;  // 0°C
# T_DROPLET = 293.15;  // 20°C
# Note: Currently only tracked, not affecting dynamics

### Hot wall (Leidenfrost regime candidate)
# T_WALL = 473.15;  // 200°C
# T_DROPLET = 293.15;  // 20°C
# Note: Full thermal coupling requires additional code

## Practical Applications

### Inkjet printing (small, fast droplets)
# D0 = 50e-6 m (50 microns)
# U0 = 5 m/s
# Grid level needs to be adjusted for length scale

### Spray coating (medium droplets)
# D0 = 100e-6 to 500e-6 m
# U0 = 1-10 m/s
# Multiple droplets require extended simulation

### Rain drop impact (large droplets)
# D0 = 5e-3 m (5 mm)
# U0 = 5-9 m/s (terminal velocity)
# Large We number, splashing expected

### Fire suppression
# D0 = 1-5 mm
# U0 = 5-15 m/s
# High momentum transfer to surface
