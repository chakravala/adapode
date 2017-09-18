! adapode, Copyright (C) 2013 Michael Reed

MODULE OdePrograms
USE Constants
USE OdeSolve
IMPLICIT NONE
CONTAINS

SUBROUTINE OdeGrid(F,filename)
IMPLICIT NONE
INTEGER::i,tol,mode
INTEGER,ALLOCATABLE::c(:)
CHARACTER::output,label
CHARACTER(LEN=*),INTENT(IN)::filename
REAL(dbl)::tmax,tmin
REAL(dbl),ALLOCATABLE::x0(:,:),xs(:),incr(:)
REAL(dbl),ALLOCATABLE::Q(:,:)
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
ALLOCATE(incr(4))
ALLOCATE(c(4))
PRINT*
PRINT*,'================================================================'
PRINT*,'Michael Reed, ODE Solver (Grid)'
PRINT*
OPEN(11,FILE=filename)
ALLOCATE(xs(4))
! Enter lower bound
READ (11,*) xs(1),xs(2),xs(3),xs(4)
! Enter upper bound
READ (11,*) incr(1),incr(2),incr(3),incr(4)
! Enter count
READ (11,*) c(1),c(2),c(3),c(4)
! Enter tmin tmax:
READ (11,*) tmin,tmax
! Enter decimal-place tolerance (e.g. 9):
READ (11,*) tol
! Enter method ID:
READ (11,*) mode
! Enter output data file name
READ (11,*) output
! Enter data file label
READ (11,*) label

 CALL Generate3dGrid(x0,xs,incr,c)
 CALL OdeSolveGrid(F,Q,x0,tmin,tmax,tol,mode)
 CALL SaveData(Q,output,filename,label)

PRINT*,'Done with calculations. Data stored in "',output,'"'
PRINT*
PRINT*,'================================================================'
PRINT*,
END SUBROUTINE OdeGrid

SUBROUTINE OdeSystem(F,filename)
IMPLICIT NONE
INTEGER::i,tol,mode
CHARACTER::output,label
CHARACTER(LEN=*),INTENT(IN)::filename
REAL(dbl)::tmax,tmin
REAL(dbl),DIMENSION(4)::x0
REAL(dbl),ALLOCATABLE::Y(:,:)
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
PRINT*
PRINT*,'================================================================'
PRINT*,'Michael Reed, ODE Solver'
PRINT*
OPEN(11,FILE=filename)
! Enter initial:
READ (11,*) x0(1),x0(2),x0(3),x0(4)
! Enter tmin tmax:
READ (11,*) tmin,tmax
! Enter decimal-place tolerance (e.g. 9):
READ (11,*) tol
! Enter method ID:
READ (11,*) mode
! Enter output data file name
READ (11,*) output
! Enter data file label
READ (11,*) label

 CALL OdeSolveSystem(F,Y,x0,tmin,tmax,tol,mode)
 CALL SaveData(Y,output,filename,label)

PRINT*,'Done with calculations. Data stored in "',output,'"'
PRINT*
PRINT*,'================================================================'
PRINT*,
END SUBROUTINE OdeSystem

SUBROUTINE ShootingSystem(F,filename)
IMPLICIT NONE
INTEGER::i,tol,mode
CHARACTER::output,label
CHARACTER(LEN=*),INTENT(IN)::filename
REAL(dbl)::tmax,tmin
REAL(dbl),DIMENSION(4)::x0
REAL(dbl),ALLOCATABLE::Y(:,:)
INTERFACE
    PURE FUNCTION F(x)
    USE Constants
    IMPLICIT NONE
    REAL(dbl),INTENT(IN),DIMENSION(:)::x
    REAL(dbl),DIMENSION(SIZE(x))::F
    END FUNCTION
END INTERFACE
PRINT*
PRINT*,'================================================================'
PRINT*,'Michael Reed, MA448, Project 1'
PRINT*
OPEN(11,FILE=filename)
! Enter initial:
READ (11,*) x0(1),x0(2),x0(3),x0(4)
! Enter tmin tmax:
READ (11,*) tmin,tmax
! Enter decimal-place tolerance (e.g. 9):
READ (11,*) tol
! Enter method ID:
READ (11,*) mode
! Enter data file output name
READ (11,*) output
! Enter data file label
READ (11,*) label

 CALL OdeSolveSystem(F,Y,x0,tmin,tmax,tol,mode)
 CALL SaveData(Y,output,filename,label)

PRINT*,'Done with calculations. Data stored in "',output,'"'
PRINT*
PRINT*,'================================================================'
PRINT*,
END SUBROUTINE ShootingSystem

END MODULE OdePrograms
