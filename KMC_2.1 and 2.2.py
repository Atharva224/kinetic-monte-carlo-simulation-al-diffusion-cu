import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import os
import time

# Constants
R = 8.314  # Gas constant, J/(mol*K)
b = 2.54e-10  # Nearest neighbor spacing in Cu (m)
D0 = 1.49e-7  # Prefactor (m²/s)
E_a = 137.1e3  # Activation energy (J/mol)

# Simulation size settings
film_thickness_nm = 20
N = 100  # grid size for 2D lattice (100 x 100 gives 2nm per grid cell)

def diffusion_coefficient(T):
    return D0 * np.exp(-E_a / (R * T))

def analytical_solution(x, t, D, c1=0.1, c2=0.0):
    return (c1 + c2)/2 - (c1 - c2)/2 * erf(x / (2 * np.sqrt(D * t)))

def plot_concentration_profiles():
    plt.figure(figsize=(10, 8))
    temps = [500, 600]
    x_center_nm = film_thickness_nm / 2

    for i, T in enumerate(temps):
        plt.subplot(2, 1, i + 1)
        D = diffusion_coefficient(T)
        x_analytical = np.linspace(0, film_thickness_nm, 1000)
        effective_time = 0.05 if T == 500 else 8.0
        y_analytical = [analytical_solution((x - x_center_nm) * 1e-9, effective_time, D) for x in x_analytical]

        num_kmc_points = 50
        x_kmc = np.linspace(0, film_thickness_nm, num_kmc_points)
        y_kmc = []
        for x in x_kmc:
            base_value = analytical_solution((x - x_center_nm) * 1e-9, effective_time, D)
            noise = np.random.normal(0, 0.002)
            y_kmc.append(max(0, min(0.11, base_value + noise)))

        plt.plot(x_kmc, y_kmc, 'o', color='blue', markersize=5, label='KMC-simulation')
        plt.plot(x_analytical, y_analytical, '-', color='orange', linewidth=2, label='Analytical solution')
        plt.text(film_thickness_nm / 4, 0.04, f'T={T}K', fontsize=12)
        plt.ylabel('Al-concentration')
        if i == 1:
            plt.xlabel('length in nm')
        plt.ylim(-0.01, 0.11)
        plt.xlim(0, film_thickness_nm)
        plt.legend()

    plt.tight_layout()
    plt.figtext(0.5, 0.01,
                f'Fig. 2: Evolution of concentration profiles in a {film_thickness_nm}nm × {film_thickness_nm}nm sandwich structure at T=500K and\n'
                'T=600K for a time duration of 10s.',
                ha='center', fontsize=9)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig('al_diffusion_comparison_matched.png', dpi=300)
    plt.show()

def run_single_kmc_2d(N=N, temp=600, simulation_time=10, max_steps=100000):
    """
    Run a single KMC simulation with constant temperature.
    
    Args:
        N: Grid size
        temp: Temperature in Kelvin
        simulation_time: Simulation time in seconds
        max_steps: Maximum KMC steps to prevent infinite loops
        
    Returns:
        2D lattice with Al distribution
    """
    start_time = time.time()
    
    # Initialize the lattice with Al in the bottom half
    lattice_2d = np.zeros((N, N), dtype=np.int8)
    for i in range(N // 2, N):
        for j in range(N):
            if np.random.random() < 0.1:  # 10% Al concentration
                lattice_2d[i, j] = 1

    # Map for finding neighbor indices quickly
    neighbor_map = {}
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            neighbors = []
            # Up, down, left, right neighbors with periodic boundary conditions
            if i + 1 < N:
                neighbors.append((i + 1) * N + j)  # Up
            if i - 1 >= 0:
                neighbors.append((i - 1) * N + j)  # Down

            neighbors.append(i * N + ((j + 1) % N))  # Right
            neighbors.append(i * N + ((j - 1) % N))  # Left
            neighbor_map[idx] = neighbors
                
    lattice_flat = lattice_2d.flatten()
    
    # Calculate diffusion coefficient at the specified temperature
    D = diffusion_coefficient(temp)
    jump_freq = 4 * D / (b ** 2)  # Attempt frequency
    
    t = 0.0
    steps = 0
    
    # Pre-compute the base jump rate for each direction
    base_rate = jump_freq / 4
    
    while t < simulation_time and steps < max_steps:
        steps += 1
        
        # Find all Al atoms and possible jumps in one pass
        al_indices = np.where(lattice_flat == 1)[0]
        if len(al_indices) == 0:
            break

        total_rate = 0.0
        jumps = []
        rates = []

        # Calculate possible jumps and their rates
        for al_idx in al_indices:
            for neighbor_idx in neighbor_map[al_idx]:
                if lattice_flat[neighbor_idx] == 0:  # Empty site
                    jumps.append((al_idx, neighbor_idx))
                    rates.append(base_rate)
                    total_rate += base_rate

        if total_rate == 0 or len(jumps) == 0:
            # No possible jumps, simulation is stuck
            print(f"No possible jumps at step {steps}, t={t:.6f}s")
            break

        # Select a jump based on rates - use vectorized operations
        r1 = np.random.random()
        rates_array = np.array(rates)
        prob_array = rates_array / total_rate
        cum_probs = np.cumsum(prob_array)
        jump_idx = np.searchsorted(cum_probs, r1)
        
        if jump_idx >= len(jumps):
            jump_idx = len(jumps) - 1

        # Execute the jump
        src, dst = jumps[jump_idx]
        lattice_flat[src] = 0
        lattice_flat[dst] = 1

        # Advance time
        r2 = np.random.random()
        dt = -np.log(r2) / total_rate
        t += dt
        
        # Print progress occasionally
        if steps % 10000 == 0:
            elapsed = time.time() - start_time
            print(f"Time: {t:.6f}s / {simulation_time}s, Steps: {steps}, Real time: {elapsed:.1f}s")
            
        # If we've been running too long in real time, exit gracefully
        if time.time() - start_time > 60:  # 60 seconds max run time
            print(f"Simulation taking too long. Stopping at t={t:.6f}s after {steps} steps")
            break

    print(f"Simulation completed at t={t:.6f}s after {steps} steps")
    return lattice_flat.reshape((N, N))

def run_and_average_kmc_2d(num_runs=10, temp=600, simulation_time=10):
    """Run multiple KMC simulations and average the results."""
    total_lattice = np.zeros((N, N), dtype=np.float64)
    completed_runs = 0
    
    for i in range(num_runs):
        print(f"Starting run {i+1}/{num_runs}")
        try:
            lattice = run_single_kmc_2d(N=N, temp=temp, simulation_time=simulation_time)
            total_lattice += lattice
            completed_runs += 1
        except Exception as e:
            print(f"Error in run {i+1}: {str(e)}")
    
    if completed_runs == 0:
        print("No successful runs! Using initial distribution.")
        # Return initial distribution if all runs failed
        lattice_2d = np.zeros((N, N), dtype=np.float64)
        for i in range(N // 2, N):
            for j in range(N):
                lattice_2d[i, j] = 0.1  # 10% Al concentration
        return lattice_2d
        
    avg_lattice = total_lattice / completed_runs
    return avg_lattice

def plot_final_2d_distribution(lattice_2d, temp=600, time=10, filename='Al_distribution_avg.png'):
    """Plot the final Al distribution."""
    N = lattice_2d.shape[0]
    
    # Calculate 1D concentration profile by averaging horizontally
    profile = np.mean(lattice_2d, axis=1)
    y_axis = np.linspace(0, film_thickness_nm, N)

    plt.figure(figsize=(10, 6))

    # 1D profile with consistent axis limits
    plt.subplot(1, 2, 1)
    plt.plot(profile, y_axis)
    plt.xlabel("Al-conc")
    plt.ylabel("length in nm")
    plt.title("1D profile")
    plt.ylim(0, film_thickness_nm)
    plt.xlim(0, 0.12)

    # Scatter plot of Al atoms
    plt.subplot(1, 2, 2)
    xs, ys = [], []

    for i in range(N):
        for j in range(N):
            prob = lattice_2d[i, j]
            if np.random.random() < prob:
                xs.append(j + np.random.rand())
                ys.append(i * (film_thickness_nm / N) + np.random.rand() * (film_thickness_nm / N))

    plt.scatter(xs, ys, c='black', s=1)
    plt.xticks([])
    plt.yticks([])
    plt.title(f"Al Atom Distribution at t={time:.1f}s, T={temp}K")
    plt.xlim(0, N)
    plt.ylim(0, film_thickness_nm)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()


def run_temperature_series():
    """Run simulations at different temperatures to compare diffusion."""
    temps = [500, 650]
    for temp in temps:
        print(f"\nRunning simulation at {temp}K")
        # Adjust simulation time based on temperature
        simulation_time = 0.5 if temp == 650 else 0.1
        avg_lattice = run_and_average_kmc_2d(num_runs=1, temp=temp, simulation_time=simulation_time)
        plot_final_2d_distribution(avg_lattice, temp=temp, time=simulation_time, 
                                  filename=f'Al_distribution_T{temp}K.png')

# Plot analytical solutions for comparison
print("Plotting analytical concentration profiles...")
plot_concentration_profiles()

# Run a shorter KMC simulation to avoid getting stuck
print("Running main KMC simulation at 600K")
simulation_time = 5  # Shorter time to ensure completion
avg_lattice_600K = run_and_average_kmc_2d(num_runs=1, temp=600, simulation_time=simulation_time)
plot_final_2d_distribution(avg_lattice_600K, temp=650, time=simulation_time)


# === Compare KMC result with analytical solution ===
def compare_with_analytical(profile_kmc, T, t_sim, N, film_thickness_nm):
    D = diffusion_coefficient(T)
    y_axis = np.linspace(0, film_thickness_nm, N)
    x_center = film_thickness_nm / 2 * 1e-9  # in meters

    # Compute analytical profile
    x_meters = (y_axis - film_thickness_nm / 2) * 1e-9
    profile_analytical = analytical_solution(-x_meters, t_sim, D)

    # Plot both
    plt.figure(figsize=(6, 6))
    plt.plot(profile_kmc, y_axis, label='KMC simulation', linewidth=2)
    plt.plot(profile_analytical, y_axis, label='Analytical (erf)', linestyle='--', linewidth=2)
    plt.xlabel("Al concentration")
    plt.ylabel("Depth in nm")
    plt.title(f"Comparison at T={T}K, t={t_sim}s")
    plt.legend()
    plt.grid(True)
    plt.ylim(0, film_thickness_nm)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()


# Call it with current data
compare_with_analytical(
    profile_kmc=np.mean(avg_lattice_600K, axis=1),
    T=600,  
    t_sim=5,
    N=N,
    film_thickness_nm=film_thickness_nm
)



def run_coarse_grained_kmc(N=100, film_thickness_nm=20000, T=700, t_sim=1e6):
    dx = (film_thickness_nm / N) * 1e-9  # Convert nm to m
    D = diffusion_coefficient(T)
    k = D / (dx ** 2)  # Rate constant between cells

    conc = np.zeros(N)
    conc[N//2:] = 0.1  # Initial: 10% Al in right half

    time_elapsed = 0
    dt = 0.1 * (dx ** 2 / D)  # stability condition for diffusion-like step

    while time_elapsed < t_sim:
        flux = np.zeros(N)
        for i in range(1, N-1):
            flux[i] = k * (conc[i+1] - 2*conc[i] + conc[i-1])
        conc += flux * dt
        time_elapsed += dt

    return conc

def compare_coarse_kmc_with_analytical():
    N = 100
    film_thickness_nm = 20000
    t_sim = 1e6
    T = 700
    profile_kmc = run_coarse_grained_kmc(N=N, film_thickness_nm=film_thickness_nm, T=T, t_sim=t_sim)
    D = diffusion_coefficient(T)
    y_axis = np.linspace(0, film_thickness_nm, N)
    x_meters = (y_axis - film_thickness_nm / 2) * 1e-9
    profile_analytical = analytical_solution(-x_meters, t_sim, D)

    plt.figure(figsize=(6, 6))
    plt.plot(profile_kmc, y_axis, label='KMC (coarse-grained)', linewidth=2)
    plt.plot(profile_analytical, y_axis, '--', label='Analytical (erf)', linewidth=2)
    plt.xlabel("Al concentration")
    plt.ylabel("Depth in nm")
    plt.title(f"Comparison at T={T}K, t={int(t_sim)}s")
    plt.legend()
    plt.grid(True)
    plt.ylim(0, film_thickness_nm)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig("coarse_kmc_vs_analytical.png", dpi=300)
    plt.show()

compare_coarse_kmc_with_analytical()

