# code developed for me 5434 - advanced intrduction to cfd
# homework 4
# ted balabanski
# dependencies
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# data reading utility
def read_matrix(filename):
    return np.genfromtxt(filename, delimiter=",", usemask=False)

# simulation details reading utility
def read_simdetails(filename):
    # open the file, read it, and split it line-by-line
    with open(filename, 'r') as f:
        raw = f.read()
    spltln = raw.split("\n")[:-1]

    # iterate through each line and create the dictionary
    res = {}
    for i in spltln:
        splt = i.split("=")
        res[splt[0]] = float(splt[1])

    return res
        

# calculate the analytical solution
def analytical(x, t, c=100, alpha=0.01, L=1):
    # calculate the two terms separately and return
    exp_term = np.exp((- alpha * (np.pi ** 2) * t) / (L ** 2))
    sin_term = np.sin(np.pi * x / L)
    return c * exp_term * sin_term
    

# declare a plotting utility
def plot_results(x, explicit, implicit, analytical, t, ax):
    # each of explicit, implicit, analytical should be 1d arrays taken at a particular time, t
    ax.plot(x, analytical, label="analytical", color='k', linestyle='--')
    ax.plot(x, explicit, label="explicit", color='r')
    ax.plot(x, implicit, label="implicit", color='b')

    # add the axes labels
    ax.set_xlabel("x [m]")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylabel("T [deg C]")
    ax.legend()
    ax.set_title(f"Time t = {t:.1f} hr")


# execute the following script when run
if __name__ == "__main__":
    # get the filenames to read as constants
    EXPLICIT_RESULTS_FILE = r"./data/explicit.csv"
    IMPLICIT_RESULTS_FILE = r"./data/implicit.csv"
    SIMDETAILS_FILE = r"./data/simdetails.txt"
    
    sim = read_simdetails(SIMDETAILS_FILE)
    explicit_results = read_matrix(EXPLICIT_RESULTS_FILE)
    implicit_results = read_matrix(IMPLICIT_RESULTS_FILE)

    # get the times to use and the corresponding indicies
    T1, T2 = 10, 100 # from the problem statement
    it1, it2 = int(T1 // sim['dt']), int(T2 // sim['dt'])
    
    # calculate the values of x from the simulation details and compute the analytical values
    x = np.linspace(0, sim['l'], explicit_results.shape[0], endpoint=True)
    analytical_t1 = analytical(x, T1, c=sim['c'], alpha=sim['alpha'], L=sim['l'])
    analytical_t2 = analytical(x, T2, c=sim['c'], alpha=sim['alpha'], L=sim['l'])

    # calculate the stability parameter
    r = (sim['alpha'] * sim['dt']) / np.power(sim['dx'], 2)
    
    # create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2))
    plot_results(x, explicit_results[:, it1], implicit_results[:, it1], analytical_t1, T1, ax1)
    plot_results(x, explicit_results[:, it2], implicit_results[:, it2], analytical_t2, T2, ax2)
    title_string = f"1D Unsteady Heat Equation Simulation\ndx = {sim['dx']:.2f} m; dt = {sim['dt']:.2f} hr; r = {r:.3f}"
    fig.suptitle(title_string)
    
    plt.savefig(f"figs/result_dx{sim['dx']*100:.0f}_dt{sim['dt']*100:.0f}.png", bbox_inches='tight', dpi=400)
    plt.show()
    
    
