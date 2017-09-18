! adapode, Copyright (C) 2013 Michael Reed

MODULE Functions
USE Constants
IMPLICIT NONE
CONTAINS


PURE FUNCTION Lorenz(x)
IMPLICIT NONE
REAL(dbl),INTENT(IN),DIMENSION(:)::x
REAL(dbl),DIMENSION(SIZE(x))::Lorenz
!  t = x(1)
! x1 = x(2)
! x2 = x(3)
! x3 = x(4)
Lorenz(1) = 1.d0
Lorenz(2) = 10.d0*(x(3) - x(2))
Lorenz(3) = x(2)*(28.d0 - x(4)) - x(3)
Lorenz(4) = x(2)*x(3) - 8.d0*x(4)/3.d0
RETURN
END FUNCTION Lorenz

PURE FUNCTION ForceFunc(rank,t)
IMPLICIT NONE
INTEGER::i
INTEGER,INTENT(IN)::rank
REAL(dbl)::x,hx
REAL(dbl),INTENT(IN)::t
REAL(dbl),DIMENSION(rank)::ForceFunc
hx = 1._dbl/rank
DO i = 1,rank
    x = i*hx
    ForceFunc(i) = 1000._dbl*(t)*DSIN(10._dbl*x + x**2)*DCOS(125._dbl*t)
END DO
END FUNCTION ForceFunc

END MODULE Functions
