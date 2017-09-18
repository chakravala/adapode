! adapode, Copyright (C) 2013 Michael Reed

MODULE OdeMultistep
USE Constants
IMPLICIT NONE
CONTAINS

PURE SUBROUTINE PredictorCorrector(F,x,xout,h,m,Fx,cp,cc,e)
USE Functions
IMPLICIT NONE
INTEGER::j,k,mp1,order
INTEGER,INTENT(INOUT)::m
REAL(dbl),INTENT(IN)::h
REAL(dbl),INTENT(OUT)::e
REAL(dbl),DIMENSION(:),INTENT(IN)::x,cp,cc
REAL(dbl),DIMENSION(SIZE(x))::s,y
REAL(dbl),DIMENSION(SIZE(x)),INTENT(OUT)::xout
REAL(dbl),DIMENSION(SIZE(cp)+1,SIZE(x)),INTENT(INOUT)::Fx
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
order = SIZE(cp)
mp1 = MOD(m,order+1) + 1
Fx(m,:) = F(x)
s = 0
DO k = 1,order
        j = MOD(m-k+order+1,order+1) + 1
        s = s + cp(k)*Fx(j,:)
END DO
y = x + h*s
Fx(mp1,:) = F(y)
s = 0
DO k = 1,order
        j = MOD(mp1-k+order+1,order+1) + 1
        s = s + cc(k)*Fx(j,:)
END DO
xout = x + h*s
e = MAXVAL(DABS( (xout-y)/xout ))
m = mp1
END SUBROUTINE PredictorCorrector

PURE SUBROUTINE AdamsBashforthMoulton2nd(m,cp,cc,Fx,eqc)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
INTEGER,INTENT(IN)::eqc
REAL(dbl),ALLOCATABLE,INTENT(OUT)::cp(:),cc(:),Fx(:,:)
m = 2
ALLOCATE(cp(m))
ALLOCATE(cc(m))
ALLOCATE(Fx(m+1,eqc))
cp = (/ 3.d0, -1.d0 /)
cc = (/ 1.d0,  1.d0 /)
cp = cp / 2.d0
cc = cc / 2.d0
END SUBROUTINE AdamsBashforthMoulton2nd

PURE SUBROUTINE AdamsBashforthMoulton3rd(m,cp,cc,Fx,eqc)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
INTEGER,INTENT(IN)::eqc
REAL(dbl),ALLOCATABLE,INTENT(OUT)::cp(:),cc(:),Fx(:,:)
m = 3
ALLOCATE(cp(m))
ALLOCATE(cc(m))
ALLOCATE(Fx(m+1,eqc))
cp = (/ 23.d0, -16.d0,  5.d0 /)
cc = (/  5.d0,   8.d0, -1.d0 /)
cp = cp / 12.d0
cc = cc / 12.d0
END SUBROUTINE AdamsBashforthMoulton3rd

PURE SUBROUTINE AdamsBashforthMoulton4th(m,cp,cc,Fx,eqc)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
INTEGER,INTENT(IN)::eqc
REAL(dbl),ALLOCATABLE,INTENT(OUT)::cp(:),cc(:),Fx(:,:)
m = 4
ALLOCATE(cp(m))
ALLOCATE(cc(m))
ALLOCATE(Fx(m+1,eqc))
cp = (/ 55.d0, -59.d0, 37.d0, -9.d0 /)
cc = (/  9.d0,  19.d0, -5.d0,  1.d0 /)
cp = cp / 24.d0
cc = cc / 24.d0
END SUBROUTINE AdamsBashforthMoulton4th

PURE SUBROUTINE AdamsBashforthMoulton5th(m,cp,cc,Fx,eqc)
IMPLICIT NONE
INTEGER,INTENT(OUT)::m
INTEGER,INTENT(IN)::eqc
REAL(dbl),ALLOCATABLE,INTENT(OUT)::cp(:),cc(:),Fx(:,:)
m = 5
ALLOCATE(cp(m))
ALLOCATE(cc(m))
ALLOCATE(Fx(m+1,eqc))
cp = (/ 1901.d0, -2774.d0, 2616.d0, -1274.d0, 251.d0 /)
cc = (/ 251.d0,    646.d0, -264.d0,   106.d0, -19.d0 /)
cp = cp / 720.d0
cc = cc / 720.d0
END SUBROUTINE AdamsBashforthMoulton5th

END MODULE OdeMultistep
