
	! Calculate x for each (e, y) pair, for each method
	DO j = 1, 9
		eps = E(j)					! eccentricity value
		DO i = 1, 32
			manom = Y(i)				! mean anomaly value	

			PRINT *,'EPSILON =',eps,', Y =',manom

			! ITERATION method
			! ------------------------------------------------------
			root = manom				! initial root
			n = 0
			CALL iterate(root, iter, tol, eps, manom, n, dec)
			! itersol(i,j) = root
			WRITE (*,7) root
7			format (1PG22.14)
			IF (n.GT.1000000) THEN
				PRINT *,'Did not converge.'
			ELSE
				PRINT *,n
			END IF

			! BISECTION method
			! ------------------------------------------------------
			lint = 0.d0					! initial left bound
			rint = 3.1415926535897932d0	! initial right bound
			n = 0
			CALL bisect(lint, rint, root, bis, tol, eps, manom, n, dec)
			! bissol(i,j) = root
			WRITE (*,8) root
8			format (1PG22.14)
			IF (n.GT.1000000) THEN
				PRINT *,'Did not converge.'
			ELSE
				PRINT *,n
			END IF

			! NEWTON method
			! ------------------------------------------------------
			root = 1.57d0			! initial root
			n = 0
			CALL iterate(root, new, tol, eps, manom, n, dec)
			! newsol(i,j) = root
			WRITE (*,9) root
9			format (1PG22.14)
			IF (n.GT.1000000) THEN
				PRINT *,'Did not converge.'
			ELSE
				PRINT *,n
			END IF

			! SECANT method
			! ------------------------------------------------------
			oldroot = 0.d0		! first initial root
			root = 1.57d0			! second initial root
			CALL iteratesecant(oldroot, root, sec, tol, eps, manom, n, dec)
			! secsol(i,j) = root
			WRITE (*,10) root
10			format (1PG22.14)
			IF (n.GT.1000000) THEN
				PRINT *,'Did not converge.'
			ELSE
				PRINT *,n
			END IF

		PRINT *,'--------------------'

		END DO
	END DO

	STOP
	END

	! ITERATION function
	! ------------------------------------------------------------
	DOUBLE PRECISION FUNCTION iter(x,eps,manom)
		implicit none
		double precision x, eps, manom
		iter = eps*dsin(x) + manom
	RETURN
	END

	! BISECTION function
	! ------------------------------------------------------------
	DOUBLE PRECISION FUNCTION bis(x, eps, manom)
		implicit none
		double precision x, eps, manom
		bis = x - eps*dsin(x) - manom
	RETURN
	END

	! NEWTON function
	! ------------------------------------------------------------
	DOUBLE PRECISION FUNCTION new(x, eps, manom)
		implicit none
		double precision x, eps, manom
		new = x + (manom-x+eps*dsin(x))/(1.0d0-eps*dcos(x))
	RETURN
	END

	! SECANT function
	! ------------------------------------------------------------
	DOUBLE PRECISION FUNCTION sec(x,z, eps, manom)
		implicit none
		double precision x, z, eps, manom
		sec = z - (z-eps*dsin(z)-manom)*(z-x)/(z-x+eps*(dsin(x)-dsin(z)))
	RETURN
	END

	! ITERATION and NEWTON subroutine
	! ------------------------------------------------------------
	SUBROUTINE iterate(root, f, tol, eps, manom, n, dec)
		implicit none
		double precision root, f, old, tol, eps, manom, new1, new2
		integer n, dec
100		old = root	
		root = f(old, eps, manom)
		IF (dec.EQ.0) THEN
			new1 = f(root, eps, manom)
			new2 = f(new1, eps, manom)
			CALL aitken(root, new1, new2, tol)
		END IF
		n = n + 1						! counts number of iterations
		IF (n.GT.1000000) GOTO 200			! stops if not converging
		IF (dabs(root-old).GE.tol) GOTO 100
200	RETURN
	END

	! BISECTION subroutine
	! ------------------------------------------------------------
	SUBROUTINE bisect(left, right, middle, f, tol, eps, manom, n, dec)
		implicit none
		double precision left, right, middle, f, fleft, fright, fmiddle, tol, eps, manom
		double precision root, old
		integer n, dec
	
		fleft = f(left, eps, manom)		! initial left bound
		fright = f(right, eps, manom)		! initial right bound

101		middle = (left+right)/2.d0			! bisect interval
		fmiddle = f(middle, eps, manom)

		IF (fleft*fmiddle.LE.0.d0) THEN	! if root to left of middle
			right = middle				! right bound becomes middle
		ELSE								! otherwise,
			left = middle					! left bound becomes middle
			fleft = fmiddle
		END IF

		old = middle
		root = (left+right)/2.d0

		n = n + 1						! counts number of iterations
		IF (n.GT.1000000) GOTO 201			! stops if not converging
		IF (dabs(left-right).GE.tol) GOTO 101 					
201	RETURN
	END

	! SECANT subroutine
	! ------------------------------------------------------------
	SUBROUTINE iteratesecant(oldroot, root, f, tol, eps, manom, n, dec)
		implicit none
		double precision oldroot, old1, root, f, tol, eps, manom, new1, new2
		integer n, dec
102		old1 = oldroot
		oldroot = root	
		root = f(old1, oldroot, eps, manom)
		IF (dec.EQ.0) THEN
			new1 = f(oldroot, root, eps, manom)
			new2 = f(root, new1, eps, manom)
			CALL aitken(root, new1, new2, tol)
		END IF
		n = n + 1						! counts number of iterations
		IF (n.GT.1000000) GOTO 202			! stops if not converging
		IF (dabs(root-oldroot).GE.tol) GOTO 102					
202	RETURN
	END

	! AITKEN'S ACCELERATION subroutine
	! ------------------------------------------------------------
	SUBROUTINE aitken(root, new1, new2, tol)
		implicit none
		double precision root, new1, new2, numer, denom, tol
		numer = root*new2 - new1*new1
		denom = new2 - 2.d0*new1 + root
		IF (denom.GT.tol) THEN
			root = numer/denom	
		END IF
203	RETURN
	END
