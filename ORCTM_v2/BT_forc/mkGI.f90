      program setup_forcingfields
      use chpar          

      real(dp),dimension(:,:),allocatable::                             &
     & UHOPE,VHOPE,TEMOUT,PRESSOUT,TDEWOUT,SWRADOUT,DWLWOUT,PRATEOUT,   &
     & U10OUT,RIV

      allocate(UHOPE(0:ie,je),VHOPE(ie,0:je),TEMOUT(ie,je),PRESSOUT(ie,je),&
     &  TDEWOUT(ie,je),SWRADOUT(ie,je),DWLWOUT(ie,je),PRATEOUT(ie,je),  &
     &  U10OUT(ie,je),RIV(ie,je) )

! forcing ----- GI*
      OPEN(52,FILE='GIWIX',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(53,FILE='GIWIY',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(54,FILE='GITEM',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(55,FILE='GITDEW', ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(56,FILE='GIPREC', ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(57,FILE='GIDWLW', ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(58,FILE='GIU10',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(59,FILE='GISWRAD',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(61,FILE='GIPRESS',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      OPEN(12,FILE='GIRIV',  ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      UHOPE=0.
      VHOPE=0.
      TEMOUT=27.1949123307819
      TDEWOUT=27.1949123307819+273.15
      PRATEOUT=0.
      U10OUT=0.
      SWRADOUT=0.
      PRESSOUT=0.
      RIV=0.

      IDAY=0
      IF3=0

      IF2=180
      WRITE(52) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int((IE+1)*JE,i8)
      WRITE(52) UHOPE
      IF2=181
      WRITE(53) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*(JE+1),i8)
      WRITE(53) VHOPE
      IF2=167
      WRITE(54) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*JE,i8)
      WRITE(54) TEMOUT
      IF2=168
      WRITE(55) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*JE,i8)
      WRITE(55) TDEWOUT
      IF2=423
      WRITE(56) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*JE,i8)
      WRITE(56) PRATEOUT
      IF2=177
      WRITE(57) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*JE,i8)
      WRITE(57) DWLWOUT
      IF2=665
      WRITE(58) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*JE,i8)
      WRITE(58) U10OUT
      IF2=276
      WRITE(59) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*JE,i8)
      WRITE(59) SWRADOUT
      IF2=151
      WRITE(61) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*JE,i8)
      WRITE(61) PRESSOUT
      IF2=999
      WRITE(12) int(IDAY,i8),int(IF2,i8),int(IF3,i8),int(IE*JE,i8)
      WRITE(12) RIV


      end program setup_forcingfields
