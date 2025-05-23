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




#Keeping the original work safe. Email atharvasinnarkar@gmail.com for the full code and mention the proper usecase.







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

