# 1D Unsteady Heat Equation Solver
Ted Balabanski

# Background

The 1D unsteady heat equation models the time-dependent temperature distribution in a homogenous solid.

$$
\begin{equation} \label{eq.2-heateqn}
  \pd{T}{t} - \alpha \pdn{T}{x}{2} = 0
\end{equation}
$$

# Discretization, BCs, and ICs

The problem is discretized by a central difference in space and backwards difference in time. The finite difference equation is given by

$$
\begin{equation} \label{eq.2-fde}
  \frac{T_i^{n+1} - T_i^n}{\Delta t} - \alpha \frac{T_{i+1}^{n+1} - 2 T_{i}^{n+1} + T_{i-1}^{n+1}}{\Delta x^2} = 0
\end{equation}
$$

The boundary conditions are $T = 0$ at both ends. The initial condition is periodic, given by

$$
\begin{equation}
	T = c \sin \left( \frac{\pi x}{L} \right)
\end{equation}
$$

# Dependencies and Use

In `main.cpp` the primary code to solve the discretized equation is presented. Two parameters may be passed, `dx`, the grid size [m] and `dt`, the time step [hr]. 

A plotting utility to plot the profile, at two time steps is provided. The analytical solution, given by the below, is also plotted

$$
	T = c \exp \left( \frac{- \alpha \pi^2 t}{L^2} \right) \sin \left( \frac{\pi x}{L} \right)
$$

To compile the binaries, simply `make` them. The code depends on the Eigen package for linear algebra routines.
