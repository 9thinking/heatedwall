#!/usr/bin/env python3
"""
Advanced analysis script for droplet impact simulation results.

This script provides:
1. Data loading and processing
2. Quantitative metrics calculation (max spreading, contact time, etc.)
3. Advanced visualization
4. Comparison with empirical correlations
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys

def load_data(filename='quantitative_output.dat'):
    """Load simulation output data."""
    try:
        data = np.loadtxt(filename, skiprows=1)
        return {
            't': data[:, 0],
            'D_D0': data[:, 1],
            'xc': data[:, 2],
            'yc': data[:, 3],
            'volume': data[:, 4],
            'KE': data[:, 5]
        }
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        print("Please run the simulation first.")
        sys.exit(1)

def calculate_metrics(data):
    """Calculate key metrics from simulation data."""
    metrics = {}

    # Maximum spreading factor
    metrics['D_max'] = np.max(data['D_D0'])
    metrics['t_max'] = data['t'][np.argmax(data['D_D0'])]

    # Contact time (when droplet center leaves contact zone)
    contact_mask = data['yc'] < 0.5  # Within half diameter of wall
    if np.any(contact_mask):
        contact_times = data['t'][contact_mask]
        if len(contact_times) > 0:
            metrics['t_contact'] = contact_times[-1] - contact_times[0]
        else:
            metrics['t_contact'] = 0.0
    else:
        metrics['t_contact'] = 0.0

    # Initial kinetic energy
    metrics['KE_initial'] = data['KE'][0] if len(data['KE']) > 0 else 0.0

    # Energy dissipation rate
    if len(data['KE']) > 10:
        metrics['energy_dissipation'] = (data['KE'][0] - data['KE'][-1]) / data['KE'][0] * 100
    else:
        metrics['energy_dissipation'] = 0.0

    # Volume conservation (should be ~constant)
    metrics['volume_change'] = (np.max(data['volume']) - np.min(data['volume'])) / np.mean(data['volume']) * 100

    return metrics

def empirical_correlation(Re, We, Oh):
    """
    Calculate maximum spreading factor using empirical correlations.

    References:
    - Pasandideh-Fard et al. (1996): D_max/D0 = sqrt((We + 12) / (3(1-cos(theta)) + 4*We/sqrt(Re)))
    - Simplified for complete wetting (theta = 0)
    """
    if Re > 0 and We > 0:
        # Pasandideh-Fard correlation (assuming complete wetting)
        D_max_empirical = np.sqrt((We + 12) / (3 + 4 * We / np.sqrt(Re)))
        return D_max_empirical
    return None

def plot_results(data, metrics, Re=None, We=None, Oh=None):
    """Create comprehensive visualization of results."""

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)

    # Plot 1: Spreading factor
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(data['t'], data['D_D0'], 'b-', linewidth=2, label='Simulation')
    ax1.axhline(y=metrics['D_max'], color='r', linestyle='--',
                label=f'Max = {metrics["D_max"]:.2f}')
    if Re and We and Oh:
        D_emp = empirical_correlation(Re, We, Oh)
        if D_emp:
            ax1.axhline(y=D_emp, color='g', linestyle=':',
                       label=f'Empirical = {D_emp:.2f}')
    ax1.set_xlabel('Time (D₀/U₀)', fontsize=11)
    ax1.set_ylabel('Spreading Factor (D/D₀)', fontsize=11)
    ax1.set_title('Droplet Spreading Evolution', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Kinetic energy
    ax2 = fig.add_subplot(gs[0, 1])
    if metrics['KE_initial'] > 0:
        KE_normalized = data['KE'] / metrics['KE_initial']
    else:
        KE_normalized = data['KE']
    ax2.plot(data['t'], KE_normalized, 'r-', linewidth=2)
    ax2.set_xlabel('Time (D₀/U₀)', fontsize=11)
    ax2.set_ylabel('Normalized Kinetic Energy', fontsize=11)
    ax2.set_title('Energy Dissipation', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # Plot 3: Center of mass trajectory
    ax3 = fig.add_subplot(gs[1, 0])
    scatter = ax3.scatter(data['xc'], data['yc'], c=data['t'],
                         cmap='viridis', s=20, alpha=0.6)
    ax3.axhline(y=0, color='k', linestyle='-', linewidth=2, label='Wall')
    ax3.set_xlabel('Radial Position (r/D₀)', fontsize=11)
    ax3.set_ylabel('Axial Position (z/D₀)', fontsize=11)
    ax3.set_title('Center of Mass Trajectory', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.set_aspect('equal', adjustable='box')
    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Time', fontsize=10)

    # Plot 4: Volume conservation
    ax4 = fig.add_subplot(gs[1, 1])
    volume_normalized = data['volume'] / np.mean(data['volume'])
    ax4.plot(data['t'], volume_normalized, 'orange', linewidth=2)
    ax4.axhline(y=1.0, color='k', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Time (D₀/U₀)', fontsize=11)
    ax4.set_ylabel('Normalized Volume', fontsize=11)
    ax4.set_title('Volume Conservation Check', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)

    # Plot 5: Phase diagram (spreading vs time)
    ax5 = fig.add_subplot(gs[2, 0])
    # Identify phases
    t_max_idx = np.argmax(data['D_D0'])
    t_spreading = data['t'][:t_max_idx+1]
    D_spreading = data['D_D0'][:t_max_idx+1]
    t_retraction = data['t'][t_max_idx:]
    D_retraction = data['D_D0'][t_max_idx:]

    ax5.plot(t_spreading, D_spreading, 'b-', linewidth=2.5, label='Spreading')
    ax5.plot(t_retraction, D_retraction, 'r-', linewidth=2.5, label='Retraction')
    ax5.axvline(x=metrics['t_max'], color='g', linestyle='--',
                label=f't_max = {metrics["t_max"]:.2f}')
    ax5.set_xlabel('Time (D₀/U₀)', fontsize=11)
    ax5.set_ylabel('D/D₀', fontsize=11)
    ax5.set_title('Spreading Phases', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    ax5.legend()

    # Plot 6: Summary metrics (text)
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.axis('off')

    summary_text = f"""
    SIMULATION SUMMARY
    {'='*40}

    Quantitative Metrics:
    • Max spreading factor: {metrics['D_max']:.3f}
    • Time to max spreading: {metrics['t_max']:.3f}
    • Contact time: {metrics['t_contact']:.3f}
    • Energy dissipation: {metrics['energy_dissipation']:.1f}%
    • Volume change: {metrics['volume_change']:.2f}%
    """

    if Re and We and Oh:
        summary_text += f"""
    Dimensionless Numbers:
    • Reynolds (Re): {Re:.1f}
    • Weber (We): {We:.1f}
    • Ohnesorge (Oh): {Oh:.4f}
    """

        # Determine regime
        if We < 5:
            regime = "Deposition"
        elif We < 30:
            regime = "Spreading"
        else:
            regime = "Splashing"
        summary_text += f"    • Impact regime: {regime}\n"

    ax6.text(0.1, 0.9, summary_text, fontsize=10, verticalalignment='top',
             fontfamily='monospace', transform=ax6.transAxes)

    plt.suptitle('Droplet Impact Analysis', fontsize=14, fontweight='bold', y=0.995)
    plt.savefig('analysis_results.png', dpi=300, bbox_inches='tight')
    print("Analysis plot saved as 'analysis_results.png'")

    return fig

def print_summary(metrics, Re=None, We=None, Oh=None):
    """Print summary to console."""
    print("\n" + "="*50)
    print("DROPLET IMPACT SIMULATION ANALYSIS")
    print("="*50)

    if Re and We and Oh:
        print(f"\nDimensionless Numbers:")
        print(f"  Reynolds number (Re): {Re:.2f}")
        print(f"  Weber number (We):    {We:.2f}")
        print(f"  Ohnesorge number (Oh): {Oh:.4f}")

        # Calculate empirical prediction
        D_emp = empirical_correlation(Re, We, Oh)
        if D_emp:
            print(f"\nEmpirical Prediction (Pasandideh-Fard):")
            print(f"  D_max/D0 = {D_emp:.3f}")

    print(f"\nSimulation Results:")
    print(f"  Maximum spreading factor: {metrics['D_max']:.3f}")
    print(f"  Time to max spreading:    {metrics['t_max']:.3f}")
    print(f"  Contact time:             {metrics['t_contact']:.3f}")
    print(f"  Energy dissipation:       {metrics['energy_dissipation']:.1f}%")
    print(f"  Volume conservation:      ±{metrics['volume_change']:.2f}%")

    if Re and We and Oh and D_emp:
        error = abs(metrics['D_max'] - D_emp) / D_emp * 100
        print(f"\nComparison with Empirical:")
        print(f"  Relative error: {error:.1f}%")

    print("\n" + "="*50 + "\n")

def main():
    """Main analysis function."""
    print("Loading simulation data...")
    data = load_data()

    print("Calculating metrics...")
    metrics = calculate_metrics(data)

    # Try to extract dimensionless numbers from log
    Re, We, Oh = None, None, None
    # These would need to be passed as arguments or read from a parameter file
    # For now, we'll calculate from default values if needed

    # Default parameters (from simulation)
    RHO_L = 1000.0
    U0 = 1.0
    D0 = 2e-3
    MU_L = 1e-3
    SIGMA = 0.072

    Re = RHO_L * U0 * D0 / MU_L
    We = RHO_L * U0**2 * D0 / SIGMA
    Oh = np.sqrt(We) / Re

    print_summary(metrics, Re, We, Oh)

    print("Generating plots...")
    plot_results(data, metrics, Re, We, Oh)

    print("Analysis complete!")

if __name__ == "__main__":
    main()
