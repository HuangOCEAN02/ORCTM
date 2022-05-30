      REAL FUNCTION SQRTRND()
      DATA B /34251./
      DATA D /0.34679201245/
      C=SQRT(B*D)
      M=INT(C)
      D=C-M
      SQRTRND=D
      END
