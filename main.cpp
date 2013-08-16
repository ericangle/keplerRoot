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
double bisect();
double secant();

double g_fpi(double x, double eps, double y);
double g_newton(double x, double eps, double y);
double g_secant(double x, double z, double eps, double y);

double f(double x, double eps, double y);
double fPrime(double x, double eps, double y);

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

/*
                        ! BISECTION method
                        ! ------------------------------------------------------
                        lint = 0.d0                                     ! initial left bound
                        rint = 3.1415926535897932d0     ! initial right bound
                        n = 0
                        CALL bisect(lint, rint, root, bis, tol, eps, manom, n, dec)
                        ! bissol(i,j) = root
                        WRITE (*,8) root
8                       format (1PG22.14)
                        IF (n.GT.1000000) THEN
                                PRINT *,'Did not converge.'
                        ELSE
                                PRINT *,n
                        END IF
*/

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


        ! SECANT subroutine
        ! ------------------------------------------------------------
        SUBROUTINE iteratesecant(oldroot, root, f, tol, eps, manom, n, dec)
                implicit none
                double precision oldroot, old1, root, f, tol, eps, manom, new1, new2
                integer n, dec
102             old1 = oldroot
                oldroot = root
                root = f(old1, oldroot, eps, manom)
                IF (dec.EQ.0) THEN
                        new1 = f(oldroot, root, eps, manom)
                        new2 = f(root, new1, eps, manom)
                        CALL aitken(root, new1, new2, tol)
                END IF
                n = n + 1                                               ! counts number of iterations
                IF (n.GT.1000000) GOTO 202                      ! stops if not converging
                IF (dabs(root-oldroot).GE.tol) GOTO 102
202     RETURN
        END


double bisect() {
  return 0;
}

/*
        ! BISECTION subroutine
        ! ------------------------------------------------------------
        SUBROUTINE bisect(left, right, middle, f, tol, eps, manom, n, dec)
                implicit none
                double precision left, right, middle, f, fleft, fright, fmiddle, tol, eps, manom
                double precision root, old
                integer n, dec

                fleft = f(left, eps, manom)             ! initial left bound
                fright = f(right, eps, manom)           ! initial right bound

101             middle = (left+right)/2.d0                      ! bisect interval
                fmiddle = f(middle, eps, manom)

                IF (fleft*fmiddle.LE.0.d0) THEN ! if root to left of middle
                        right = middle                          ! right bound becomes middle
                ELSE                                                            ! otherwise,
                        left = middle                                   ! left bound becomes middle
                        fleft = fmiddle
                END IF

                old = middle
                root = (left+right)/2.d0

                n = n + 1                                               ! counts number of iterations
                IF (n.GT.1000000) GOTO 201                      ! stops if not converging
                IF (dabs(left-right).GE.tol) GOTO 101
201     RETURN
        END
*/

double secant() {
  return 0;
}

double g_fpi(double x, double eps, double y) {
  return y + eps*sin(x);
}

double g_newton(double x, double eps, double y) {
  return x + (y - x + eps*sin(x))/(1.0 - eps*cos(x));
}

double f(double x, double eps, double y) {
  return x - eps*sin(x) - y;
}

double fPrime(double x, double eps, double y) {
  return 1.0 - eps*cos(x);
}

double g_secant(double x_n, double x_nmm, double eps, double y) {
  return x_n - (x_n - x_nmm)*f(x_n,eps,y)/(f(x_n,eps,y) - f(x_nmm,eps,y));
}

/* Aitken's acceleration
SUBROUTINE aitken(root, new1, new2, tol)
  implicit none
  double precision root, new1, new2, numer, denom, tol
  numer = root*new2 - new1*new1
  denom = new2 - 2.d0*new1 + root
  IF (denom.GT.tol) THEN
    root = numer/denom
  END IF
*/
