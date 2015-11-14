#include <iostream>
#include <cmath>
using namespace std;

// Can probably combine iterate and iterate secant
double iterate(double g(double x, double eps, double y), double root, double eps, double y, double tol, int& n, int stop);
double iterateSecant(double g(double x_n, double x_nmm, double eps, double y), double root, double oldRoot, double eps, double y, double tol, int& n, int stop);
double bisect(double f(double x, double eps, double y), double left, double right, double eps, double y, double tol, int& n, int stop);

double g_fpi(double x, double eps, double y);
double g_newton(double x, double eps, double y);
double g_secant(double x, double z, double eps, double y);

double f(double x, double eps, double y);
double fPrime(double x, double eps, double y);

void outputLine(string algType, double root, int n, int stop);

int main() {
  double pi = 3.14159265358979323846;
  double piOverTwo = pi/2.0;  

  // Tolerance: want roots accurate to 14 decimal places
  double tol = 0.5e-14;

  // Eccentricity values
  double E [] = {0.0, 0.1, 0.2, 0.4, 0.8, 0.9, 0.95, 0.975, 0.9875};
 
  double eps, y, root;
  int n;
  int stop = 1000000; // stop after this many steps

  for (int i = 0; i < 9; i++) {
    eps = E[i];
    for (int j = 0; j < 32; j++) {
      y = ((double) j)/10.0; // Mean anomaly values: 0.0, 0.1, ..., 3.1
      cout << "eps = " << eps << ", y = " << y << endl;
      cout << "----------------------" << endl;      

      // Fixed point iteration
      n = 0;
      root = iterate(g_fpi, piOverTwo, eps, y, tol, n, stop);
      outputLine("FPI", root, n, stop);

      // Newton
      n = 0;
      root = iterate(g_newton, piOverTwo, eps, y, tol, n, stop);
      outputLine("NEW", root, n, stop);

      // Interval bisection
      n = 0;
      root = bisect(f, 0.0, pi, eps, y, tol, n, stop);
      outputLine("BIS", root, n, stop);

      // Secant
      n = 0;
      root = iterateSecant(g_secant, piOverTwo, 0.0, eps, y, tol, n, stop);
      outputLine("SEC", root, n, stop);      

      cout << endl;
    }    
  }

  return 0;
}
  
double iterate(double g(double x, double eps, double y), double root, double eps, double y, double tol, int& n, int stop) {
  double error = tol + 1.0;
  double old;
  while (error > tol && n < stop) {
    old = root;
    root = g(old, eps, y);
    error = abs(root - old);
    n++;
  }
  return root;
}

double iterateSecant(double g(double x_n, double x_nmm, double eps, double y), double root, double oldRoot, double eps, double y, double tol, int& n, int stop) {
  double error = tol + 1.0;
  double old, older;
  while (error > tol && n < stop) {
    old = root;
    older = oldRoot;
    root = g(old, older, eps, y);
    oldRoot = old;
    error = abs(root - old);
    n++;
  }
  return root;
}

double bisect(double f(double x, double eps, double y), double left, double right, double eps, double y, double tol, int& n, int stop) {
  double error = tol + 1.0;
  double root;
  double fleft, fright, middle, fmiddle;

  while (error > tol && n < stop) {
    fleft = f(left, eps, y);      // initial left bound
    fright = f(right, eps, y);    // initial right bound

    middle = (left + right)/2.0;      // bisect interval
    fmiddle = f(middle, eps, y);

    if (fleft*fmiddle <= 0.0) {   // if root to left of middle
      right = middle;              // right bound becomes middle
    }
    else {                                              
      left = middle;               // left bound becomes middle
    }
    root = (left + right)/2.0;
    error = abs(left - right);
    n++;
  }
  return root;
}

double g_fpi(double x, double eps, double y) {
  return x - f(x, eps, y);
}

double g_newton(double x, double eps, double y) {
  return x - f(x, eps, y)/fPrime(x, eps, y);
}

double g_secant(double x_n, double x_nmm, double eps, double y) {
  return x_n - (x_n - x_nmm)*f(x_n,eps,y)/(f(x_n,eps,y) - f(x_nmm,eps,y));
}

double f(double x, double eps, double y) {
  return x - eps*sin(x) - y;
}

double fPrime(double x, double eps, double y) {
  return 1.0 - eps*cos(x);
}

void outputLine(string algType, double root, int n, int stop) {
  if (n == stop) {
    cout << algType << ": Did not converge.";
  }
  else {
    cout << algType << ": " << root << '\t' << n;
  }
  cout << endl;
}
