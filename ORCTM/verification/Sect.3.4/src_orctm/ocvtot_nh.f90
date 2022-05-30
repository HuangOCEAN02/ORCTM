      SUBROUTINE OCVTOT_NH
! Update the 2D Poisson problem for vertical bottom boundary condition -2021.06.12 by H.H in Qingdao
#ifdef NON_HYDROSTATIC
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      use mo_petsc
      use mo_obcs
            
! calculate GU,GV,GW
      DO 3171 K=1,KE
      DO 3171 J=1,JE
      DO 3171 I=I_start,IE
        GU(I,J,K)=(UOO(I,J,K)-UOO_1(I,J,K))/DT
 3171  CONTINUE

      DO 3172 K=1,KE
      DO 3172 J=J_start,JE
      DO 3172 I=1,IE
        GV(I,J,K)=(VOE(I,J,K)-VOE_1(I,J,K))/DT
 3172  CONTINUE

      DO 3173 K=1,KEP
      DO 3173 J=1,JE
      DO 3173 I=1,IE
        GW(I,J,K)=(WNO(I,J,K)-WNO_1(I,J,K))/DT
 3173  CONTINUE

! calculate del(G) -->> divergence of (GU,GV,GW)
      DO 3174 K=1,KE
      DO 3174 J=1,JE
      DO 3174 I=(I_start+1),IE
         DIVGU(I,J,K)=(GU(I,J,K)-GU(I-1,J,K))/DLXP(I,J)
3174  CONTINUE
      CALL bounds_exch(DIVGU)

      DO 3175 K=1,KE
      DO 3175 J=(J_start+1),JE
      DO 3175 I=1,IE
         DIVGV(I,J,K)=(GV(I,J-1,K)-GV(I,J,K))/DLYP(I,J)

3175  CONTINUE
      CALL bounds_exch(DIVGV)

      DO 3176 K=2,KE
      DO 3176 J=1,JE
      DO 3176 I=1,IE
         DIVGW(I,J,K)=(GW(I,J,K)-GW(I,J,K+1))/max(DDPO(I,J,K),1.E-3)

3176  CONTINUE

      DO 3177 J=1,JE
      DO 3177 I=1,IE
         DIVGW(I,J,1)=(GW(I,J,1)-GW(I,J,2))/(max(DDPO(I,J,1),1.E-3)+ZO(I,J))
3177  CONTINUE
      CALL bounds_exch(DIVGW)

      DO 3170 K=1,KE
      DO 3170 J=1,JE
      DO 3170 I=1,IE
        DIVG(I,J,K)=(DIVGU(I,J,K)+DIVGV(I,J,K)+DIVGW(I,J,K))
3170  CONTINUE

#ifdef CG2D4NONHY  
! Extract the GU,GV,GW at bottom   
      DO 3200 K=1,KE
      DO 3200 J=1,JE
      DO 3200 I=I_start,IE
        if (K .eq. 1) then
          if(AMSUO(I,J,K) .eq. 0)  BOTGU(I,J) = 0
        elseif (K .eq. KE) then
            if (AMSUO(I,J,K) .ne. 0) then
               BOTGU(I,J)=GU(I,J,K)
            endif
        else
           if (AMSUO(I,J,K-1) .ne. 0 .AND. AMSUO(I,J,K) .eq. 0 ) then
               BOTGU(I,J)=GU(I,J,K-1)
           endif
         endif
3200  CONTINUE
      

      DO 3201 K=1,KE
      DO 3201 J=J_start,JE
      DO 3201 I=1,IE
       if (K .eq. 1) then
         if(AMSUE(I,J,K) .eq. 0)  BOTGV(I,J) = 0
       elseif (K .eq. KE) then
           if (AMSUE(I,J,K) .ne. 0) then
            BOTGV(I,J)=GV(I,J,K)
           endif
       else
           if (AMSUE(I,J,K-1) .ne. 0 .AND. AMSUE(I,J,K) .eq. 0 ) then
               BOTGV(I,J)=GV(I,J,K-1)
           endif
         endif
3201  CONTINUE
 
! Attention! The vertical velocity at bottom is zero!!! 
      DO 3202 K=1,KE
      DO 3202 J=1,JE
      DO 3202 I=1,IE
      
        if (K .eq. 1) then
          if(WETO(I,J,K) .eq. 0)  BOTGW(I,J) = 0
        elseif (K .eq. KE) then
           if (WETO(I,J,K) .ne. 0) then
            BOTGW(I,J)=GW(I,J,K)/DDPO(I,J,K)
           endif
        else
         if (WETO(I,J,K-1) .ne. 0 .AND. WETO(I,J,K) .eq. 0 ) then
               BOTGW(I,J)=GW(I,J,K-1)/DDPO(I,J,K-1)
         endif
         endif              
 3202  CONTINUE
 
 
! calculate del(G) -->> divergence of (GU,GV,GW) at bottom
      DO 3203 J=1,JE
      DO 3203 I=(I_start+1),IE
         DIVBOTGU(I,J)=(BOTGU(I,J)-BOTGU(I-1,J))/DLXP(I,J)
3203  CONTINUE
      CALL bounds_exch(DIVBOTGU)

      DO 3204 J=(J_start+1),JE
      DO 3204 I=1,IE
         DIVBOTGV(I,J)=(BOTGV(I,J-1)-BOTGV(I,J))/DLYP(I,J)
3204  CONTINUE
      CALL bounds_exch(DIVBOTGV)

      DO 3205 J=1,JE
      DO 3205 I=1,IE   
         DIVBOTGW(I,J)=BOTGW(I,J)
3205  CONTINUE 
      CALL bounds_exch(DIVBOTGW)
 
      DO 3206 J=1,JE
      DO 3206 I=1,IE     
         DIVG_bot_surf(I,J) = (DIVBOTGU(I,J)+DIVBOTGV(I,J)+DIVBOTGW(I,J))   &
        &  - DIVG(I,J,1)*weto(I,J,1)
3206  CONTINUE      
  
      
! Solve 2D Poisson Problem 
      call petsc_solve_2D
      CALL bounds_exch(PNH_bot)
#endif /*CG2D4NONHY*/      
            
! Solve 3D Poisson Problem
      call petsc_solve_3D
      CALL bounds_exch(PNH)
      
! Open Boundary conditions for nonhydrostatic pressure
      if (ICYCLI_X.EQ.0) then
        if (have_g_is .and. west) then
          DO J=2,JE1
          DO K=1,KE
            PNH(1,J,K)=PNH(2,J,K)
          ENDDO
          ENDDO
        endif
        if (have_g_ie .and. east) then
          DO J=2,JE1
          DO K=1,KE
            PNH(IE,J,K)=PNH(IE1,J,K)
          ENDDO
          ENDDO
        endif
      endif

      if (ICYCLI_Y.EQ.0) then
        if (have_g_js .and. north) then
          DO I=1,IE
          DO K=1,KE
            PNH(I,1,K)=PNH(I,2,K)
          ENDDO
          ENDDO
        endif
        if (have_g_je .and. south) then
          DO I=1,IE
          DO K=1,KE
            PNH(I,JE,K)=PNH(I,JE1,K)
          ENDDO
          ENDDO
        endif
      endif

! correct u,v with pnh
      DO 1995 K=1,KE
      DO 1995 J=2,JE1
      DO 1995 I=(I_start+1),IE1
        UOO(I,J,K)=UOO(I,J,K)+AMSUO(I,J,K)*DTDXUO(I,J)                  &
                                          *(PNH(I,J,K)-PNH(I+1,J,K))
1995  CONTINUE
      CALL bounds_exch(UOO)

      DO 1996 K=1,KE
      DO 1996 J=(J_start+1),JE1
      DO 1996 I=2,IE1
        VOE(I,J,K)=VOE(I,J,K)+AMSUE(I,J,K)*DPYE(I,J)                    &
                                          *(PNH(I,J+1,K)-PNH(I,J,K))
1996  CONTINUE
      CALL bounds_exch(VOE)
!
!RJ       CALL CONTRO(221)
!
!======================================================================
!
!     B)
!
!     VERTICAL VELOCITY = VERTICAL INTEGRAL OF DIVERGENCE OF
!                             HORIZONTAL VELOCITY FIELD
!

      DO 7811 K=KEP,1,-1
      DO 7811 J=1,JE
      DO 7811 I=1,IE
      WNO(I,J,K) = ZERO
 7811 CONTINUE

      DO 781 K=KE,1,-1
      DO 781 J=(J_start+1),JE
      DO 781 I=(I_start+1),IE

      WNO(I,J,K) = WNO(I,J,K+1)                                         &
     & + DTI * WETO(I,J,K) * (                                          &
     & DTDXPO(I,J)   * (   UOO(I-1,J,K) * DDUO(I-1,J,K)*DLYU(I-1,J)     &
     &  - UOO(I,J,K)   * DDUO(I,J,K)*DLYU(I,J)   )/DLYP(I,J)  )         &
     & + DTI*WETO(I,J,K)* (                                             &
     & + DTDYO(I,J) *                                                   &
     & (   VOE(I,J,K)     * DDUE(I,J,K)*DLXV(I,J)                       &
     & - VOE(I,J-1,K)   * DDUE(I,J-1,K)*DLXV(I,J-1)  )/DLXP(I,J)  )

781    CONTINUE

      CALL bounds_exch(WNO)

      DO 17 K=1,KE
      DO 17 J=1,JE
      DO 17 I=I_start,IE
       UOO_1(I,J,K)=UOO(I,J,K)
17    CONTINUE

      DO 18 K=1,KE
      DO 18 J=J_start,JE
      DO 18 I=1,IE
       VOE_1(I,J,K)=VOE(I,J,K)
18    CONTINUE

      DO 19 K=1,KEP
      DO 19 J=1,JE
      DO 19 I=1,IE
       WNO_1(I,J,K)=WNO(I,J,K)
       WO(I,J,K)=WNO(I,J,K)
19    CONTINUE

#endif
      RETURN
      ENDSUBROUTINE OCVTOT_NH

