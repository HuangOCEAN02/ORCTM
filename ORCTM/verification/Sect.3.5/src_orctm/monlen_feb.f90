      SUBROUTINE MONLEN_FEB(NACYEAR,LENGTH)
            LENGTH = 28
            IF (    (MOD(NACYEAR,  4).EQ.0 .AND.                        &
     &          MOD(NACYEAR,100).NE.0)                                  &
     &          .OR. MOD(NACYEAR,400).EQ.0  ) LENGTH = 29
      RETURN
      END
