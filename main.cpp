// The Kepler Problem and Root Finding

// g++ main.cpp -o executable_name

// Note about PDF explaining background

// This program solves Kepler's equation y = x - e*sin(x) for the eccentric anomaly
// x, given the mean anomaly y and eccentricity e, using an iteration method and the
// bisection, Newton, and secant root-finding methods.

#include <iostream>
#include <math.h>
#include <cmath>
using namespace std;


double iterate(double g(double x, double eps, double y), double root, double eps, double y, double tol, int& n);
double newton();
double bisect();
double secant();

double g_fpi(double x, double eps, double y);
double g_newton(double x, double eps, double y);

int main() {
  // Tolerance: want roots accurate to 14 decimal places
  double tol = 0.5e-14;

  // Eccentricity values
  double E [] = {0.0, 0.1, 0.2, 0.4, 0.8, 0.9, 0.95, 0.975, 0.9875};

  // Mean anomaly values: 0.0, 0.1, ..., 3.1
  double Y [32];
  for (int i = 0; i < 32; i++) {
    Y[i] = ((double) i)/10.0;
  }
 
  double eps, y, root;
  int n;

  for (int i = 0; i < 9; i++) {
    eps = E[i];
    for (int j = 0; j < 32; j++) {
      y = Y[j];

      cout << "eps = " << eps << ", y = " << y << endl;

      // Fixed point iteration
      n = 0;
      root = y;
      root = iterate(g_fpi, root, eps, y, tol, n);
      if (n == 1000000) {
        cout << "Did not converge." << endl;
      }
      else {
        cout << root << '\t' << n << endl; 
      }

      // Newton
      n = 0;
      root = iterate(g_fpi, root, eps, y, tol, n);
      if (n == 1000000) {
        cout << "Did not converge." << endl;
      }
      else {
        cout << root << '\t' << n << endl;
      }

      // Interval bisection
      n = 0;
      root = iterate(g_fpi, y, eps, y, tol, n);
      if (n == 1000000) {
        cout << "Did not converge." << endl;
      }
      else {
        cout << root << '\t' << n << endl;
      }

      // Secant
      n = 0;
      root = iterate(g_fpi, y, eps, y, tol, n);
      if (n == 1000000) {
        cout << "Did not converge." << endl;
      }
      else {
        cout << root << '\t' << n << endl;
      }

    }    
  }

  return 0;
}
  
double iterate(double g(double x, double eps, double y), double root, double eps, double y, double tol, int& n) {
  double error = tol + 1.0;
  double old;
  while (error > tol && n < 1000000) {
    old = root;
    root = g(old, eps, y);
    error = abs(root - old);
    n++;
  }
  return root;
}

double g_fpi(double x, double eps, double y) {
  return y + eps*sin(x);
}

double g_newton(double x, double eps, double y) {
  return x + (y - x + eps*sin(x))/(1.0 - eps*cos(x));
}

double newton() {
  return 0;
}

double bisect() {
  return 0;
}

double secant() {
  return 0;
}
