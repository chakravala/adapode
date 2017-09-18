! adapode, Copyright (C) 2013 Michael Reed

MODULE Constants
IMPLICIT NONE
INTEGER,PARAMETER::dbl = KIND(1.d0)
REAL(dbl),PARAMETER::PI = 16._dbl*DATAN(0.2_dbl) + &
                            4._dbl*DATAN(1._dbl/239._dbl)
CONTAINS

SUBROUTINE SaveData(Y,filename,header,label)
IMPLICIT NONE
INTEGER::i,j
CHARACTER(LEN=*),INTENT(IN)::filename,header,label
REAL(dbl),DIMENSION(:,:),INTENT(IN)::Y
OPEN(12,FILE=filename)
WRITE(12,*) '# Michael Reed, ',header
WRITE(12,*) '# ',label
DO i=1,SIZE(Y,1)
        WRITE(12,*) (Y(i,j), j = 1,SIZE(Y,2))
END DO
CLOSE(12)
END SUBROUTINE SaveData
END MODULE Constants
