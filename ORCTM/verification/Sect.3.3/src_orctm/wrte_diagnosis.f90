      SUBROUTINE WRTE_DIAGNOSIS(KDTDAY,KDAYS,KMONTS,KYEARS,KMEAN,NANF)
#ifdef DIAGNOSIS
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_DIAGNOSIS

      
!   *KYEARS*   - actual year.
!   *KMONT*    - actual month.
!   *KDAYS*    - actual day.
      
      NDTDAY=NINT(86400./DT)
      
!spy 5-day average
      IF (KMEAN.EQ.5) THEN
        NSTEP=0
        DO I=1,(KMONTS-1)
          NSTEP=NSTEP+MONLEN(I)*NDTDAY
        ENDDO
        NSTEP=NSTEP+(KDAYS-1)*NDTDAY+KDTDAY
        NEND=NDTDAY*5
        NANF=MOD((NSTEP-1),NEND)+1

!HH step by step
      ELSEIF (KMEAN.EQ.4) THEN

        NANF=1
        NEND=1
      
!HH Daily average
      ELSEIF (KMEAN.EQ.1) THEN

        NANF=KDTDAY
        NEND=NDTDAY
      
!HH Monthly average
      ELSEIF (KMEAN.EQ.2) THEN

        IF (KDAYS.EQ.1) THEN
          NANF=KDTDAY
        ELSEIF (KDAYS.GT.1) THEN
          NANF=KDTDAY+((KDAYS-1)*NDTDAY)
        ENDIF
        NEND=NDTDAY*MONLEN(KMONTS)
 
!HH Yearly average
      ELSEIF (KMEAN.EQ.3) THEN

        IF ((KDAYS.EQ.1).AND.(KMONTS.EQ.1)) THEN
          NANF=KDTDAY
        ELSE
          NANF=NANF+1
        ENDIF
        NEND=0
        DO I=1,12
          NEND=NEND+MONLEN(I)
        ENDDO
        NEND=NEND*NDTDAY

!HH No diagnostic output
      ELSEIF (KMEAN.EQ.0) THEN
        CONTINUE
      ELSE
         STOP 'Stop in WRTE_DIAGNOSIS due to wrong parameter.'
      ENDIF
 
      IF (KMEAN.NE.0) THEN
       WRITE(IO_STDOUT,*) 'Diag module nanf,nend',nanf,nend  
            
 ! spy: IVEC: hrizontal grid -1 -> w point;  0 -> p point; 1 -> u point; 2 -> v point   
! 3D temporal means
         CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Tt &
     &              ,THO,Tt_point,0)    
         CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Ss &
     &              ,SAO,Ss_point,0)    
         CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Uu &
     &              ,UKO,Uu_point,1)    
         CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Vv &
     &              ,VKE,Vv_point,2)
      
#ifdef NON_HYDROSTATIC

        CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Ww  &
     &              ,WNO,Ww_point,-1)    
        CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Pn  &
     &              ,PNH,Pn_point,-1)           
!        CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Div &
!     &              ,DIVG,div_point,-1)

#else

         CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Ww &
     &              ,WO,Ww_point,-1)
      
#endif

        CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Po &
     &              ,PO,Po_point,0)      
 
#ifdef DIFFDIAG

        CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Av  &
     &              ,AVO,Av_point,-1) 
        CALL MMEAN_DIAG3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Dv  &
     &              ,DVO,Dv_point,-1)
  
#endif

! 2D temporal means
         CALL MMEAN_DIAG2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_Zz &
     &              ,ZO,Zz_point)  
#ifdef BAROTROPIC
         CALL MMEAN_DIAG2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_US &
     &              ,USO,Us_point)  
         CALL MMEAN_DIAG2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_DIAG_VS &
     &              ,VSE,Vz_point)       
#endif
     
      ENDIF 
#endif
     RETURN    
    END SUBROUTINE WRTE_DIAGNOSIS    
