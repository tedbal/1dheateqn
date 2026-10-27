// code developed for me 5434 - advanced intrduction to cfd
// homework 4
// ted balabanski
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <numeric>
#include <cstdlib>
#include <array>
#include <bits/stdc++.h>
#include <numbers>
#include <fstream>

// include the eigen package for linear algebra routines
#include <Eigen/Dense>

// initialize and declare a default for the time length of the simulation
// declare defaults for the timestep and grid spacing
constexpr double TMAX { 150.0 }; // hours
constexpr double L { 1.0 }; // m
double dt { 0.1 }; // hours
double dx { 0.1 }; // m

// declare constants for the simulation
double TA { 0.0 }; // boundary condition A, deg C
double TB { 0.0 }; // boundary condition B, deg C
double ALPHA { 0.01 }; // thermal diffusivity, m2/hr
double C { 100.0 }; // initial condition coeff, deg C

// function: intialize the temperature array for results
Eigen::MatrixXd initializeTempArray(int nx, int nt, float dx) {
  // declare and calculate the initial temperature
  Eigen::MatrixXd T (nx, nt);
  for (int i { 1 }; i < nx - 1; ++i) {
    T(i, 0) = C * std::sin(M_PI * dx * i / L);
  }

  // set boundary conditions
  T(0, Eigen::indexing::all) = Eigen::MatrixXd::Constant(1, nt, TA);
  T(nx - 1, Eigen::indexing::all) = Eigen::MatrixXd::Constant(1, nt, TB);

  return T;
}

// save out the simulation details
void saveSim(const std::string& filename, float dx, float dt) {
  // save the file
  std::ofstream file(filename);

  // write out parameter values
  file << "tmax=" << TMAX << "\n";
  file << "l=" << L << "\n";
  file << "dx=" << dx << "\n";
  file << "dt=" << dt << "\n";
  file << "c=" << C << "\n";
  file << "alpha=" << ALPHA << "\n";

  file.close();
}

// save out the results to a file (csv format)
void saveMatrix(const std::string& filename, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& matrix) {
  // open the file
  std::ofstream file(filename);
  
  // write out the matrix as a csv format
  for (Eigen::Index i = 0; i < matrix.rows(); ++i) {
    for (Eigen::Index j = 0; j < matrix.cols(); ++j) {
      if (j > 0) file << ',';
      file << matrix(i, j);	
      }
      file << '\n';
    }

  file.close();
}
  

int main(int argc, char *argv[]) {
  // allow the user to pass in a time step and grid spacing 
  if (argc > 2) {
    dx = std::stod( argv[1] );
    dt = std::stod( argv[2] );
  }

  // calculate the stability parameter, r
  double r = (ALPHA * dt) / std::pow(dx, 2);
  
  // calculate the number of temporal points to resolve and create the mesh
  int NX = static_cast<int>(std::round(L / dx));
  int NT = static_cast<int>(std::round(TMAX / dt));
  
  // initialize two arrays, one for the explicit integration and one for the implicit integration
  Eigen::MatrixXd T_explicit = initializeTempArray(NX, NT, dx);
  Eigen::MatrixXd T_implicit = initializeTempArray(NX, NT, dx);
  std::cout << "T rows: " << T_explicit.rows() << ", cols: " << T_explicit.cols() << std::endl;

  // EXPLICIT METHOD
  // iterate through the temporal dimension (marching) [omit n = 0 b/c ICs have been set]
  for (int n { 0 }; n < NT - 1; ++n) {
    // iterate through the spatial dimension [omit i = 0, NX-1 b/c BCs have been set]
    for (int i { 1 }; i < NX - 1; ++i) {
      T_explicit(i, n+1) = r * T_explicit(i-1, n) - (1 - 2*r) * T_explicit(i, n) + r * T_explicit(i+1, n);
    }
  }
  

  // IMPLICIT METHOD  
  // the inverse problem to be solved for the implicit time integration is A T^{n+1} = T^n at time step n
  // since the grid spacing and timestep do not change, the matrix A is constant and can be declared outside the loop
  Eigen::MatrixXd A (NX, NX);
  std::cout << "A rows: " << A.rows() << ", cols: " << A.cols() << std::endl;
  
  // assign the boundary conditions, then (off-)diagonal terms
  A(0, 0) = 1.0;
  A(NX-1, NX-1) = 1.0;
  for (int i { 1 }; i < NX - 1; ++i) {
    A(i, i) = -2.0 * r - 1.0; // diagonal term
    A(i, i-1) = r; // lower off-diagonal term
    A(i, i+1) = r; // upper off-diagonal term
  }

  // begin the integration in time
  for (int n { 0 }; n < NT-1; ++n) {
    // solve the inverse problem using SVD decomposition and set the solution to the next timestep
    T_implicit(Eigen::indexing::all, n+1) = A.bdcSvd(Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV).solve(-1 * T_implicit(Eigen::indexing::all, n));
  }
  
  // SAVE OUT THE RESULTS
  saveSim("./data/simdetails.txt", dx, dt);
  saveMatrix("./data/explicit.csv", T_explicit);
  saveMatrix("./data/implicit.csv", T_implicit);

  return 0;
}
