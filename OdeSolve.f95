! adapode, Copyright (C) 2013 Michael Reed

MODULE OdeSolve
USE Constants
USE Functions
IMPLICIT NONE
CONTAINS

SUBROUTINE OdeSolveGrid(F,Q,x0,a,b,tol,mode)
IMPLICIT NONE
INTEGER::s,v,j
INTEGER,INTENT(IN)::tol,mode
REAL(dbl),INTENT(IN)::a,b
REAL(dbl),DIMENSION(:,:),INTENT(IN)::x0
REAL(dbl),ALLOCATABLE::Y(:,:)
REAL(dbl),ALLOCATABLE,INTENT(OUT)::Q(:,:)
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
ALLOCATE(Q(0,SIZE(x0,2)))
DO j = 1,SIZE(x0,1)
        CALL OdeSolveSystem(F,Y,x0(j,:),a,b,tol,mode)
        s = SIZE(Q,1)
        v = SIZE(Y,1)
        CALL ResizeArray(Q,s,v,SIZE(x0,2))
        Q(s+1:s+v,:) = Y(:,:)
END DO
END SUBROUTINE OdeSolveGrid

SUBROUTINE OdeSolveSystem(F,Y,x0,a,b,tol,mode)
USE OdeEuler
USE OdeRungeKutta
USE OdeMultistep
IMPLICIT NONE
INTEGER::i,j,m,method,mm,iflag,q,itmax,initialize,N
INTEGER,INTENT(IN)::tol,mode
REAL(dbl)::h,d,hmin,hmax,emin,emax,t,e
REAL(dbl),INTENT(IN)::a,b
REAL(dbl),DIMENSION(:),INTENT(IN)::x0
REAL(dbl),DIMENSION(SIZE(x0))::x,xout
REAL(dbl),ALLOCATABLE,INTENT(OUT)::Y(:,:)
REAL(dbl),ALLOCATABLE::c(:),brs(:,:),ce(:),cp(:),cc(:),Y2(:,:),Fx(:,:)
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
h = 2.d0**(-tol)
N = (b-a)*2.d0**tol
emax = 10.d0**(-tol)
emin = 10.d0**(-tol-3)
hmin = 1.d-16
hmax = 1.d-4
itmax = 5**(tol+5)
IF (ALLOCATED(Y) .eqv. .TRUE.) THEN
DEALLOCATE(Y)
END IF
ALLOCATE(Y(N+1,SIZE(x0)))
Y(1,:) = x0 ! Set initial conditions
IF (mode == 1) THEN
        ! Explicit Euler Method
        DO i = 1,N
                x = Y(i,:)
                Y(i+1,:) = ExplEuler(F,x,h)
        END DO
        method = 0
ELSEIF (mode == 2) THEN
        ! Implicit Euler Method
        DO i = 1,N
                x = Y(i,:)
                Y(i+1,:) = ImplEuler(F,x,h)
        END DO
        method = 0
ELSEIF (mode == 3) THEN
        ! Improved Euler, Heun's Method
        DO i = 1,N
                x = Y(i,:)
                Y(i+1,:) = ImprovedHeun(F,x,h)
        END DO
        method = 0
ELSEIF (mode == 4) THEN
        ! Midpoint 2nd Order RK Method
        CALL Midpoint2nd(m,brs,c)
        method = 1
ELSEIF (mode == 5) THEN
        ! Kutta's 3rd Order RK Method
        CALL Kutta3rd(m,brs,c)
        method = 1
ELSEIF (mode == 6) THEN
        ! Classical 4th Order RK Method
        CALL Classical4th(m,brs,c)
        method = 1
ELSEIF (mode == 11) THEN
        ! Adaptive Heun-Euler
        CALL AdaptiveHeunEuler(m,brs,c,ce)
        method = 2
ELSEIF (mode == 12) THEN
        ! Adaptive Bogacki-Shampine RK23 Method
        CALL BogackiShampineRK23(m,brs,c,ce)
        method = 2
ELSEIF (mode == 13) THEN
        ! Adaptive Fehlberg RK45 Method
        CALL FehlbergRK45(m,brs,c,ce)
        method = 2
ELSEIF (mode == 14) THEN
        ! Adaptive Cash-Karp RK45 Method
        CALL CashKarpRK45(m,brs,c,ce)
        method = 2
ELSEIF (mode == 15) THEN
        ! Adaptive Dormand-Prince RK45 Method
        CALL DormandPrinceRK45(m,brs,c,ce)
        method = 2
ELSEIF (mode == 22) THEN
        ! Adams-Bashorth-Moulton 2nd Order
        CALL AdamsBashforthMoulton2nd(mm,cp,cc,Fx,SIZE(x0))
        ! Classical RK4 Initializer
        CALL Classical4th(m,brs,c)
        method = 3
ELSEIF (mode == 23) THEN
        ! Adams-Bashorth-Moulton 3rd Order
        CALL AdamsBashforthMoulton3rd(mm,cp,cc,Fx,SIZE(x0))
        ! Classical RK4 Initializer
        CALL Classical4th(m,brs,c)
        method = 3
ELSEIF (mode == 24) THEN
        ! Adams-Bashorth-Moulton 4th Order
        CALL AdamsBashforthMoulton4th(mm,cp,cc,Fx,SIZE(x0))
        ! Classical RK4 Initializer
        CALL Classical4th(m,brs,c)
        method = 3
ELSEIF (mode == 25) THEN
        ! Adams-Bashorth-Moulton 5th Order
        CALL AdamsBashforthMoulton5th(mm,cp,cc,Fx,SIZE(x0))
        ! Classical RK4 Initializer
        CALL Classical4th(m,brs,c)
        method = 3
ELSEIF (mode == 32) THEN
        ! Adaptive Adams-Bashorth-Moulton 2nd Order
        CALL AdamsBashforthMoulton2nd(mm,cp,cc,Fx,SIZE(x0))
        ! Classical RK4 Initializer
        CALL Classical4th(m,brs,c)
        method = 4
ELSEIF (mode == 33) THEN
        ! Adaptive Adams-Bashorth-Moulton 3rd Order
        CALL AdamsBashforthMoulton3rd(mm,cp,cc,Fx,SIZE(x0))
        ! Classical RK4 Initializer
        CALL Classical4th(m,brs,c)
        method = 4
ELSEIF (mode == 34) THEN
        ! Adaptive Adams-Bashorth-Moulton 4th Order
        CALL AdamsBashforthMoulton4th(mm,cp,cc,Fx,SIZE(x0))
        ! Classical RK4 Initializer
        CALL Classical4th(m,brs,c)
        method = 4
ELSEIF (mode == 35) THEN
        ! Adaptive Adams-Bashorth-Moulton 5th Order
        CALL AdamsBashforthMoulton5th(mm,cp,cc,Fx,SIZE(x0))
        ! Classical RK4 Initializer
        CALL Classical4th(m,brs,c)
        method = 4
END IF
IF (method == 1) THEN ! Single Step
        DO i = 1,N
                x = Y(i,:)
                Y(i+1,:) = GeneralMethod(F,x,h,m,brs,c)
        END DO
ELSEIF (method == 2) THEN ! Adaptive Single Step
        i = 0
        iflag = 1
        CALL CheckhSize(h,hmin,hmax)
        DO WHILE (i .LE. itmax)
                i = i + 1
                x = Y(i,:)
                CALL CheckFinalStep(x(1),b,h,iflag)
                CALL ResizeArray(Y,i,10000,SIZE(x0))
                CALL AdaptiveMethod(F,x,xout,h,m,brs,c,ce,e)
                Y(i+1,:) = xout
                IF (iflag == 0) THEN
                        EXIT
                END IF
                IF (e < emin) THEN
                        h = 2*h
                END IF
                IF (e > emax) THEN
                        h = h/2.d0
                        i = i - 1
                END IF
                CALL CheckhSize(h,hmin,hmax)
                CALL ShowProgress(x,i,b)
        END DO
        CALL TruncateArray(Y,i,SIZE(x0))
ELSEIF (method == 3) THEN ! Multistep
        DO i = 1,(mm-1)
                x = Y(i,:)
                Y(i+1,:) = GeneralMethod(F,x,h,m,brs,c)
                Fx(i,:) = F(x)
        END DO
        q = mm
        DO i = mm,N
                x = Y(i,:)
                CALL PredictorCorrector(F,x,xout,h,q,Fx,cp,cc,e)
                Y(i+1,:) = xout
        END DO
ELSEIF (method == 4) THEN ! Adaptive Multistep
        i = 0
        iflag = 1
        initialize = 0
        CALL CheckhSize(h,hmin,hmax)
        DO WHILE (i .LE. itmax)
                i = i + 1
                x = Y(i,:)
                CALL CheckFinalStep(x(1),b,h,iflag)
                CALL ResizeArray(Y,i+mm,10000,SIZE(x0))
                IF (initialize == 0) THEN
                        DO j = 1,mm-1
                                Y(i+j,:) = GeneralMethod(F,x,h,m,brs,c)
                                Fx(j,:) = F(x)
                                x = Y(i+j,:)
                        END DO
                        q = mm
                        initialize = 1
                        i = i + mm - 1
                END IF
                CALL PredictorCorrector(F,x,xout,h,q,Fx,cp,cc,e)
                Y(i+1,:) = xout
                IF (iflag == 0) THEN
                        EXIT
                END IF
                IF (e < emin) THEN
                        h = 2*h
                        i = i - CEILING(mm/2.d0)
                        initialize = 0
                END IF
                IF (e > emax) THEN
                        h = h/2.d0
                        i = i - CEILING(mm/2.d0)
                        initialize = 0
                END IF
                CALL CheckhSize(h,hmin,hmax)
                !CALL ShowProgress(x,i,b)
        END DO
        CALL TruncateArray(Y,i,SIZE(x0))
END IF
END SUBROUTINE OdeSolveSystem

PURE SUBROUTINE CheckhSize(h,hmin,hmax)
IMPLICIT NONE
REAL(dbl),INTENT(IN)::hmin,hmax
REAL(dbl),INTENT(INOUT)::h
IF (DABS(h) < hmin) THEN
        h = DSIGN(hmin,h)
END IF
IF (DABS(h) > hmax) THEN
        h = DSIGN(hmax,h)
END IF
END SUBROUTINE CheckhSize

PURE SUBROUTINE ResizeArray(Y,i,increment,eqc)
IMPLICIT NONE
INTEGER,INTENT(IN)::i,increment,eqc
REAL(dbl),ALLOCATABLE::Y2(:,:)
REAL(dbl),ALLOCATABLE,INTENT(INOUT)::Y(:,:)
IF (SIZE(Y,1) < i+1) THEN
        ALLOCATE(Y2(i+increment,eqc))
        Y2(1:i,:) = Y(1:i,:)
        DEALLOCATE(Y)
        ALLOCATE(Y(i+increment,eqc))
        Y(1:i,:) = Y2(1:i,:)
        DEALLOCATE(Y2)
END IF
END SUBROUTINE ResizeArray

PURE SUBROUTINE TruncateArray(Y,i,eqc)
IMPLICIT NONE
INTEGER,INTENT(IN)::i,eqc
REAL(dbl),ALLOCATABLE::Y2(:,:)
REAL(dbl),ALLOCATABLE,INTENT(INOUT)::Y(:,:)
IF (SIZE(Y,1) > i+1) THEN
        ALLOCATE(Y2(i,eqc))
        Y2(1:i,:) = Y(1:i,:)
        DEALLOCATE(Y)
        ALLOCATE(Y(i,eqc))
        Y = Y2
        DEALLOCATE(Y2)
END IF
END SUBROUTINE TruncateArray

PURE SUBROUTINE CheckFinalStep(a,b,h,iflag)
IMPLICIT NONE
INTEGER,INTENT(OUT)::iflag
REAL(dbl)::d
REAL(dbl),INTENT(IN)::a,b
REAL(dbl),INTENT(INOUT)::h
d = DABS(b-a)
IF (d .LE. DABS(h)) THEN
        iflag = 0
        h = DSIGN(d,h)
END IF
END SUBROUTINE CheckFinalStep

SUBROUTINE ShowProgress(x,i,b)
IMPLICIT NONE
INTEGER,INTENT(IN)::i
REAL(dbl),INTENT(IN)::b
REAL(dbl),DIMENSION(:),INTENT(IN)::x
IF (Mod(i,75000) == 11) THEN
        PRINT*,x(1)," out of ",b
END IF
END SUBROUTINE ShowProgress

PURE SUBROUTINE Generate3dGrid(x0,xL,xU,c)
IMPLICIT NONE
INTEGER::k,j,l,n,loc
INTEGER,DIMENSION(4),INTENT(IN)::c
REAL(dbl),DIMENSION(4)::incr
REAL(dbl),DIMENSION(4),INTENT(INOUT)::xL,xU
REAL(dbl),ALLOCATABLE,INTENT(OUT)::x0(:,:)
k = 1
DO j = 1,4
        IF (c(j) .GE. 2) THEN
                incr(j) = (xU(j) - xL(j)) / (c(j)-1)
        ELSE
                incr(j) = 0.d0
        END IF
        k = k*c(j)
END DO
ALLOCATE(x0(k,4))
DO j = 1,c(2)
        DO l = 1,c(3)
                DO n = 1,c(4)
                        loc = n + ((c(3))*(l-1)) + (((c(4))**2)*(j-1))
                        x0(loc,1) = xL(1)
                        x0(loc,2) = xL(2) + incr(2)*(j-1)
                        x0(loc,3) = xL(3) + incr(3)*(l-1)
                        x0(loc,4) = xL(4) + incr(4)*(n-1)
                        !PRINT*,x0(loc,1),x0(loc,2),x0(loc,3),x0(loc,4)
                END DO
        END DO
END DO
END SUBROUTINE Generate3dGrid

END MODULE OdeSolve
