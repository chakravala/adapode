! adapode, Copyright (C) 2013 Michael Reed

MODULE OdeRungeKutta
USE Constants
USE Functions
IMPLICIT NONE
CONTAINS

PURE SUBROUTINE AdaptiveMethod(F,x,xout,h,m,b,c,ce,e)
IMPLICIT NONE
INTEGER::r
INTEGER,INTENT(IN)::m
REAL(dbl),INTENT(IN)::h
REAL(dbl),INTENT(OUT)::e
REAL(dbl),DIMENSION(:),INTENT(IN)::x,c,ce
REAL(dbl),DIMENSION(:,:),INTENT(IN)::b
REAL(dbl),DIMENSION(SIZE(x))::p,s
REAL(dbl),DIMENSION(SIZE(x)),INTENT(OUT)::xout
REAL(dbl),DIMENSION(m,SIZE(x))::K
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
K = ExplicitK(F,x,h,b,m)
p = c(1)*K(1,:)
s = ce(1)*K(1,:)
DO r = 2,m
        p = p + c(r)*K(r,:)
        s = s + ce(r)*K(r,:)
END DO
e = MAXVAL(DABS(h*s))
xout = x + h*p
END SUBROUTINE AdaptiveMethod

PURE FUNCTION GeneralMethod(F,x,h,m,b,c)
IMPLICIT NONE
INTEGER::r
INTEGER,INTENT(IN)::m
REAL(dbl),INTENT(IN)::h
REAL(dbl),DIMENSION(:),INTENT(IN)::x,c
REAL(dbl),DIMENSION(:,:),INTENT(IN)::b
REAL(dbl),DIMENSION(SIZE(x))::p,GeneralMethod
REAL(dbl),DIMENSION(m,SIZE(x))::K
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
K = ExplicitK(F,x,h,b,m)
p = c(1)*K(1,:)
DO r = 2,m
        p = p + c(r)*K(r,:)
END DO
GeneralMethod = x + h*p
END FUNCTION GeneralMethod

PURE FUNCTION ExplicitK(F,x,h,b,m)
USE Functions
IMPLICIT NONE
INTEGER::r,s
INTEGER,INTENT(IN)::m
REAL(dbl),INTENT(IN)::h
REAL(dbl),DIMENSION(:,:),INTENT(IN)::b
REAL(dbl),DIMENSIOn(:),INTENT(IN)::x
REAL(dbl),DIMENSION(SIZE(x))::z
REAL(dbl),DIMENSION(m,SIZE(x))::ExplicitK
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
ExplicitK(1,:) = F(x)
DO r = 2,m
        z = 0.d0
        DO s = 1,r-1
                z = z + b(r,s)*ExplicitK(s,:)
        END DO
        ExplicitK(r,:) = F(x + h*z)
END DO
RETURN
END FUNCTION ExplicitK

PURE SUBROUTINE Midpoint2nd(m,brs,c)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
REAL(dbl),ALLOCATABLE,INTENT(OUT)::c(:),brs(:,:)
m = 2
ALLOCATE(c(m))
ALLOCATE(brs(m,m))
c = (/ 0.d0, 1.d0 /)
brs(1,:) = (/ 0.d0 ,0.d0 /)
brs(2,:) = (/ 0.5d0,0.d0 /)
END SUBROUTINE Midpoint2nd

PURE SUBROUTINE Kutta3rd(m,brs,c)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
REAL(dbl),ALLOCATABLE,INTENT(OUT)::c(:),brs(:,:)
m = 3
ALLOCATE(c(m))
ALLOCATE(brs(m,m))
c = (/ 1.d0, 4.d0, 1.d0/)
c = c / 6.d0
brs(1,:) = (/ 0.d0 , 0.d0, 0.d0 /)
brs(2,:) = (/ 0.5d0, 0.d0, 0.d0 /)
brs(3,:) = (/ -1.d0, 2.d0, 0.d0 /)
END SUBROUTINE Kutta3rd

PURE SUBROUTINE Classical4th(m,brs,c)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
REAL(dbl),ALLOCATABLE,INTENT(OUT)::c(:),brs(:,:)
m = 4
ALLOCATE(c(m))
ALLOCATE(brs(m,m))
c = (/ 1.d0, 2.d0, 2.d0, 1.d0 /)
c = c / 6.d0
brs(1,:) = (/ 0.d0 ,0.d0 ,0.d0,0.d0 /)
brs(2,:) = (/ 0.5d0,0.d0 ,0.d0,0.d0 /)
brs(3,:) = (/ 0.d0 ,0.5d0,0.d0,0.d0 /)
brs(4,:) = (/ 0.d0 ,0.d0 ,1.d0,0.d0 /)
END SUBROUTINE Classical4th

PURE SUBROUTINE AdaptiveHeunEuler(m,brs,c,ce)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
REAL(dbl),ALLOCATABLE,INTENT(OUT)::c(:),ce(:),brs(:,:)
m = 2
ALLOCATE(c(m))
ALLOCATE(ce(m))
ALLOCATE(brs(m,m))
ce = (/ 0.5d0 ,0.5d0 /)
c =  (/ 1.d0  ,0.d0  /)
ce = c-ce ! convert constants to compute error term
brs(1,:) = (/ 0.d0 ,0.d0 /)
brs(2,:) = (/ 1.d0 ,0.d0 /)
END SUBROUTINE AdaptiveHeunEuler

PURE SUBROUTINE BogackiShampineRK23(m,brs,c,ce)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
REAL(dbl),ALLOCATABLE,INTENT(OUT)::c(:),ce(:),brs(:,:)
m = 4
ALLOCATE(c(m))
ALLOCATE(ce(m))
ALLOCATE(brs(m,m))
ce = (/ 2.d0/9.d0,  1.d0/3.d0, 4.d0/9.d0, 0.d0      /)
c =  (/ 7.d0/24.d0, 1.d0/4.d0, 1.d0/3.d0, 1.d0/8.d0 /)
ce = c-ce ! convert constants to compute error term
brs(1,:) = (/ 0.d0      ,0.d0      ,0.d0      ,0.d0 /)
brs(2,:) = (/ 0.5d0     ,0.d0      ,0.d0      ,0.d0 /)
brs(3,:) = (/ 0.d0      ,0.75d0    ,0.d0      ,0.d0 /)
brs(4,:) = (/ 2.d0/9.d0 ,1.d0/3.d0 ,4.d0/9.d0 ,0.d0/)
END SUBROUTINE BogackiShampineRK23

PURE SUBROUTINE FehlbergRK45(m,brs,c,ce)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
REAL(dbl),ALLOCATABLE,INTENT(OUT)::c(:),ce(:),brs(:,:)
m = 6
ALLOCATE(c(m))
ALLOCATE(ce(m))
ALLOCATE(brs(m,m))
ce = (/ 25.d0/216.d0, 0.d0, 1408.d0/2565.d0, 2197.d0/4104.d0, -0.2d0, 0.d0 /)
c = (/ 16.d0/135.d0, 0.d0, 6656.d0/12825.d0, 28561.d0/56430.d0, &
        -0.18d0, 2.d0/55.d0 /)
ce = c-ce ! convert constants to compute error term
brs(1,:) = (/ 0.d0            ,0.d0             ,0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(2,:) = (/ 0.25d0          ,0.d0             ,0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(3,:) = (/ 0.09375d0       ,0.28125d0        ,0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(4,:) = (/ 1932.d0/2197.d0 ,-7200.d0/2197.d0 ,7296.d0/2197.d0 ,&
              0.d0            ,0.d0             ,0.d0 /)
brs(5,:) = (/ 439.d0/216.d0   ,-8.d0            ,3680.d0/513.d0  ,&
              -845.d0/4104.d0 ,0.d0             ,0.d0/)
brs(6,:) = (/ -8.d0/27.d0     ,2.d0             ,3544.d0/2565.d0 ,&
              1859.d0/4104.d0 ,-11.d0/40.d0     ,0.d0/)
END SUBROUTINE FehlbergRK45

PURE SUBROUTINE CashKarpRK45(m,brs,c,ce)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
REAL(dbl),ALLOCATABLE,INTENT(OUT)::c(:),ce(:),brs(:,:)
m = 6
ALLOCATE(c(m))
ALLOCATE(ce(m))
ALLOCATE(brs(m,m))
ce = (/ 37.d0/378.d0, 0.d0, 250.d0/621.d0, 125.d0/594.d0, 0d0, &
        512.d0/1771.d0 /)
c = (/ 2825.d0/27648.d0, 0.d0, 18575.d0/48384.d0, 13525.d0/55296.d0, &
        277.d0/14336.d0, 0.25d0 /)
ce = c-ce ! convert constants to compute error term
brs(1,:) = (/ 0.d0              ,0.d0           ,0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(2,:) = (/ 0.2d0             ,0.d0           ,0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(3,:) = (/ 3.d0/40.d0        ,9.d0/40.d0     ,0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(4,:) = (/ 3.d0/40.d0        ,-9.d0/10.d0    ,6.d0/5.d0  ,&
              0.d0              ,0.d0           ,0.d0/)
brs(5,:) = (/ -11.d0/54.d0      ,5.d0/2.d0      ,-70.d0/27.d0 ,&
              35.d0/27.d0       ,0.d0           ,0.d0/)
brs(6,:) = (/ 1631.d0/55296.d0  ,175.d0/512.d0  ,575.d0/13824.d0 ,&
              44275.d0/110592.d0,253.d0/4096.d0 ,0.d0/)
END SUBROUTINE CashKarpRK45

PURE SUBROUTINE DormandPrinceRK45(m,brs,c,ce)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
REAL(dbl),ALLOCATABLE,INTENT(OUT)::c(:),ce(:),brs(:,:)
m = 7
ALLOCATE(c(m))
ALLOCATE(ce(m))
ALLOCATE(brs(m,m))
ce = (/ 5179.d0/57600.d0, 0.d0, 7571.d0/16695.d0, 393.d0/640.d0, &
        -92097.d0/339200.d0, 187.d0/2100.d0, 1.d0/40.d0 /)
c = (/ 35.d0/384.d0, 0.d0, 500.d0/1113.d0, 125.d0/192.d0, &
        -2187.d0/6784.d0, 11.d0/84.d0, 0.d0 /)
ce = c-ce ! convert constants to compute error term
brs(1,:) = (/ 0.d0            ,0.d0      , 0.d0, 0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(2,:) = (/ 0.2d0           ,0.d0      , 0.d0, 0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(3,:) = (/ 3.d0/40.d0      ,9.d0/40.d0, 0.d0, 0.d0 ,0.d0 ,0.d0, 0.d0 /)
brs(4,:) = (/ 44.d0/45.d0     ,-56.d0/15.d0     ,32.d0/9.d0 ,&
              0.d0            ,0.d0             ,0.d0 ,0.d0/)
brs(5,:) = (/ 19372.d0/6561.d0,-25360.d0/2187.d0,64448.d0/6561.d0 ,&
              -212.d0/729.d0  ,0.d0             ,0.d0 ,0.d0/)
brs(6,:) = (/ 9017.d0/3168.d0 ,-355.d0/33.d0    ,46732.d0/5247.d0 ,&
              49.d0/176.d0    ,-5103.d0/18656.d0,0.d0 ,0.d0/)
brs(7,:) = (/ 35.d0/384.d0    ,0.d0             ,500.d0/1113.d0 ,&
              125.d0/192.d0   ,-2187.d0/6784.d0 ,11.d0/84.d0,0.d0 /)
END SUBROUTINE DormandPrinceRK45

END MODULE OdeRungeKutta
