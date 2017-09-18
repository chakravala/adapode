! adapode, Copyright (C) 2013 Michael Reed

PROGRAM compute
USE Constants
USE OdePrograms
IMPLICIT NONE
INTEGER::i
REAL(dbl),DIMENSION(4)::b,x,a,d,c

PRINT*
PRINT*,"1. ODE Solver (Lorenz System)"
PRINT*,"2. ODE Grid   (Lorenz System)"
PRINT*
PRINT*,"Enter an integer:"
READ*,i
PRINT*

IF (i == 1) THEN
 CALL OdeSystem(Lorenz,'LorenzSystem')
ELSE IF (i == 2) THEN
 CALL OdeGrid(Lorenz,'LorenzGrid')
END IF

END PROGRAM compute
