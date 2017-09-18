! adapode, Copyright (C) 2013 Michael Reed

MODULE OdeEuler
USE Functions
USE Constants
IMPLICIT NONE
CONTAINS

PURE FUNCTION ExplEuler(F,x,h)
IMPLICIT NONE
REAL(dbl),INTENT(IN)::h
REAL(dbl),DIMENSION(:),INTENT(IN)::x
REAL(dbl),DIMENSION(SIZE(x))::ExplEuler
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
ExplEuler = x + h*F(x)
END FUNCTION ExplEuler

PURE FUNCTION ImplEuler(F,x,h)
IMPLICIT NONE
REAL(dbl),INTENT(IN)::h
REAL(dbl),DIMENSION(:),INTENT(IN)::x
REAL(dbl),DIMENSION(SIZE(x))::ImplEuler
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
! to be implemented
ImplEuler = x + h*F( x + h*F(x) )
END FUNCTION ImplEuler

PURE FUNCTION ImprovedHeun(F,x,h)
IMPLICIT NONE
REAL(dbl),INTENT(IN)::h
REAL(dbl),DIMENSION(:),INTENT(IN)::x
REAL(dbl),DIMENSION(SIZE(x))::ImprovedHeun,Fx
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
Fx = F(x)
ImprovedHeun = x + 0.5d0*h*( Fx + F( x + h*Fx ) )
END FUNCTION ImprovedHeun

END MODULE OdeEuler
