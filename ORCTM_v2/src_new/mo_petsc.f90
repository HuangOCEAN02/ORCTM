MODULE MO_PETSC
! Solve Poisson Equation for non-hydrostatic assumption by Krylov subspace method
! Compute the matrix and right-hand-side vector that define the linear system, Ax = b. by Spy
! Update the Open Boundary Conditions for Solving Poission Problem     -2021.04.28 by H.H in Qingdao
! Update the 2D Poisson problem for vertical bottom boundary condition -2021.06.12 by H.H in Qingdao
! Update the Periodic Boundary Condtitions                             -2021.06.28 by H.H in Qingdao
#ifdef NON_HYDROSTATIC

#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscsys.h>

! Include "petscksp.h" so that we can use KSP solvers.  Note that this file
!  automatically includes:
!     petscsys.h       - base PETSc routines   petscvec.h - vectors
!     petscmat.h - matrices
!     petscis.h     - index sets            petscksp.h - Krylov subspace methods
!     petscviewer.h - viewers               petscpc.h  - preconditioners  
  use petscdmda
  use petscksp
  use petscdm
  use petscsys

  USE MO_PARAM1
  USE MO_COMMO1
  USE MO_PARALLEL
  USE MO_OBCS
! Declare the Vec rhs(right-hand-side) & res(solution) for Petsc usage
  DM   da
  KSP  ksp
  Vec rhs,res,NullVecs
  
#ifdef CG2D4NONHY  
  DM   da_2d
  KSP  ksp_2d
  Vec rhs_2d,res_2d
#endif /*CG2D4NONHY*/

  PetscInt   ::izero = 0
  PetscInt   ::ione = 1
  PetscScalar::zero_scalar = 0.0 
  PetscBool  ::Neumann = .FALSE.
  
  CONTAINS

    SUBROUTINE petsc_start
 
      PetscInt ie_g2,je_g2      
      PetscInt,pointer::limit_x(:),limit_y(:),size_x(:),size_y(:),size_z(:)
      allocate(limit_x(0:nprocx),limit_y(0:nprocy))
      allocate(size_x(nprocx),size_y(nprocy),size_z(1))

      ie_g2=ie_g-2
      je_g2=je_g-2
      size_z(1)=ke
! Processors allocation
      DO i=0,nprocx
        limit_x(i)=2+i*ie_g2/nprocx
        if(i.ge.1) then
          size_x(i)=limit_x(i)-limit_x(i-1)
        endif
      ENDDO
      DO j=0,nprocy
        limit_y(j)=2+j*je_g2/nprocy
        if(j.ge.1) then
          size_y(j)=limit_y(j)-limit_y(j-1)
        endif
      ENDDO

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      if (ierr .ne. 0) then
        print*,'Unable to initialize PETSc'
        stop
      endif
      
! Create the 3D-KSP&DMDA
      call KSPCreate(MPI_COMM_WORLD,ksp,ierr)
      call KSPSetOptionsPrefix(ksp,'cg3d_', ierr)

#ifdef CYCLICWE           
       call DMDACreate3d(PETSC_COMM_WORLD,&   
       & DM_BOUNDARY_PERIODIC,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,&        
       & DMDA_STENCIL_BOX,ie_g2,je_g2,ke, &
       & nprocx,nprocy,1,1,1,size_x,size_y,size_z,da,ierr)  
#endif 
#ifdef CYCLICNS
       call DMDACreate3d(PETSC_COMM_WORLD,&   
       & DM_BOUNDARY_NONE,DM_BOUNDARY_PERIODIC,DM_BOUNDARY_NONE,& 
       & DMDA_STENCIL_BOX,ie_g2,je_g2,ke, &
       & nprocx,nprocy,1,1,1,size_x,size_y,size_z,da,ierr)     
#endif 

#ifndef CYCLICWE 
#ifndef CYCLICNS
       call DMDACreate3d(PETSC_COMM_WORLD,& 
       & DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,&
       & DMDA_STENCIL_BOX,ie_g2,je_g2,ke, &
       & nprocx,nprocy,1,1,1,size_x,size_y,size_z,da,ierr)                
#endif     
#endif    
      call DMSetFromOptions(da,ierr)      
      call DMSetUp(da,ierr)
      call DMDASetInterpolationType(da,DMDA_Q0,ierr)   
      
      call KSPSetDM(ksp,da,ierr)   
      call KSPSetComputeOperators(ksp,compute_matrix,iopt,ierr)
      call KSPSetFromOptions(ksp,ierr)
      call KSPSetUp(ksp,ierr)
      call DMCreateGlobalVector(da,rhs,ierr)
      call VecDuplicate(rhs,res,ierr)
      
#ifdef CG2D4NONHY  
! Create the 2D-KSP&DMDA
      call KSPCreate(MPI_COMM_WORLD,ksp_2d,ierr)     
      call KSPSetOptionsPrefix(ksp_2d,'cg2d_', ierr)
      call DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,&
           DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,ie_g2,je_g2, &
           nprocx,nprocy,1,1,size_x,size_y,da_2d,ierr)          
      call DMSetOptionsPrefix(da_2d,'cg2d_',ierr)
      call DMSetFromOptions(da_2d,ierr)        
      call DMSetUp(da_2d,ierr)    
      call DMDASetInterpolationType(da_2d,DMDA_Q0,ierr) 
      
      call KSPSetDM(ksp_2d,da_2d,ierr)      
      call KSPSetComputeOperators(ksp_2d,compute_matrix_2d,iopt,ierr)   
      call KSPSetFromOptions(ksp_2d,ierr)
      call KSPSetUp(ksp_2d,ierr)                 
      call DMCreateGlobalVector(da_2d,rhs_2d,ierr)
      call VecDuplicate(rhs_2d,res_2d,ierr)     

#endif /*CG2D4NONHY*/
        
    END SUBROUTINE petsc_start


    SUBROUTINE petsc_end
    
      call KSPDestroy(ksp,ierr)
      call DMDestroy(da,ierr)    
      
#ifdef CG2D4NONHY      
      call KSPDestroy(ksp_2d,ierr)
      call DMDestroy(da_2d,ierr)
#endif /*CG2D4NONHY*/          

      call PetscFinalize(ierr)
      
    END SUBROUTINE petsc_end



    SUBROUTINE petsc_solve_2D
#ifdef CG2D4NONHY      
      PetscInt   xs,ys,xm,ym      
      PetscScalar,pointer::res_2d_array(:,:)     
      
      call DMDAGetCorners(da_2d,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
      call KSPSetComputeRHS(ksp_2d,ComputeRHS_2d,ctx,ierr);
      call KSPSolve(ksp_2d,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr) 
      call KSPGetSolution(ksp_2d,res_2d,ierr)  
      call DMDAVecGetArrayReadF90(da_2d,res_2d,res_2d_array,ierr)   
      do 30 j=2,je1
      do 30 i=2,ie1
        PNH_bot(i,j)=res_2d_array(xs+i-2,ys+j-2)
30    continue

      call DMDAVecRestoreArrayReadF90(da_2d,res_2d,res_2d_array,ierr) 
      call KSPSetInitialGuessNonzero(ksp_2d,PETSC_TRUE,ierr)       
#endif /*CG2D4NONHY*/               
    END SUBROUTINE petsc_solve_2D  
    
    SUBROUTINE ComputeRHS_2d(ksp_2d,bb,ctx,ierr) 
#ifdef CG2D4NONHY        
      DM           da_2d
      KSP          ksp_2d
      Vec          bb
      PetscInt     ctx,n
      MatNullSpace nullspace
      PetscInt   xs,ys,xm,ym
      PetscErrorCode  ierr
      PetscScalar,pointer::rhs_2d_array(:,:)
      PetscScalar,pointer::Nullvecs_array(:)
      

      call KSPGetDM(ksp_2d,da_2d,ierr)
      call DMDAGetCorners(da_2d,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
      print*, xs,ys,xm,ym
      call DMDAVecGetArrayF90(da_2d,bb,rhs_2d_array,ierr)
      do 20 j=2,je1
      do 20 i=2,ie1
        rhs_2d_array(xs+i-2,ys+j-2)=DIVG_bot_surf(i,j)
20    continue  

      call DMDAVecRestoreArrayF90(da_2d,bb,rhs_2d_array,ierr)
       
      call VecAssemblyBegin(bb,ierr)
      call VecAssemblyEnd(bb,ierr)
     
       if (Neumann)  then

         call VecCreateSeq(PETSC_COMM_SELF,ione,NullVecs,ierr)
         call VecGetArrayF90(NullVecs,NullVecs_array,ierr)
            NullVecs_array(ione)=zero_scalar
         call VecRestoreArrayF90(NullVecs,NullVecs_array,ierr)      
         call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,izero,NullVecs,nullspace,ierr)
         call MatNullSpaceRemove(nullspace,bb,NullVecs,ierr)
         call MatNullSpaceDestroy(nullspace,ierr)
         
       endif
          
     return
#endif /*CG2D4NONHY*/      
    END SUBROUTINE ComputeRHS_2d


    SUBROUTINE petsc_solve_3D
      PetscInt   xs,ys,zs,xm,ym,zm,its,maxits
      PetscReal  norm,rtol,abstol,dtol
      PetscScalar,pointer::res_array(:,:,:)     
      
      call DMDAGetCorners(da,xs,ys,zs,xm,ym,zm,ierr)
      call KSPSetComputeRHS(ksp,ComputeRHS,ctx,ierr)   
! Prescribe the convergence tolerance
! the defaults are rtol=1e-5, atol=1e-50, dtol=1e+5, and maxits=1e+4.
     ! call KSPSetTolerances(ksp,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION, &
     !      PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)  
           
!      call KSPSetConvergenceTest(ksp,MyKSPConverged, &
!           PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)           
!      call KSPGetTolerances(ksp,rtol,abstol,dtol,maxits)  
!      if (p_pe==p_io) then
!            write(IO_STDOUT,*) '  rtol:',rtol
!            write(IO_STDOUT,*) 'abstol:',abtol
!            write(IO_STDOUT,*) '  dtol:',dtol
!            write(IO_STDOUT,*) 'maxits:',maxits
!      endif         
  
      call KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr) 
      call KSPGetSolution(ksp,res,ierr)    
!      call KSPSetConvergenceTest(ksp,MyKSPConverged, &
!           PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)  
      call DMDAVecGetArrayReadF90(da,res,res_array,ierr)   
      do 30 k=1,ke
      do 30 j=2,je1
      do 30 i=2,ie1
        PNH(i,j,k)=res_array(xs+i-2,ys+j-2,zs+k-1)
        
!        if (K .GE. 104 .and. K .LE. 111 ) then
!          if (abs(PNH(i,j,k)) .GE. (0.8 - (K-104)*0.1)) then
!                PNH(i,j,k)=(0.8 - (K - 104)*0.1)*PNH(i,j,K)/abs(PNH(i,j,k))
!          endif
!        endif

!        if (K .GT. 111) then
!          if (abs(PNH(i,j,k)) .GE. 0.1) then
!                PNH(i,j,k)=0.1*PNH(i,j,K)/abs(PNH(i,j,k))
!          endif
!        endif

!       if (K .GE. 104) then
!          PNH(I,J,K) = PNH(I,J,103)
!       endif
      
30    continue

      call DMDAVecRestoreArrayReadF90(da,res,res_array,ierr)
! Tell the iterative solver that the initial guess is nonzero
! Set previous solution as initial guess for next solve.    
      call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)    
    END SUBROUTINE petsc_solve_3D  



    SUBROUTINE ComputeRHS(ksp,bb,ctx,ierr) 
      DM           da
      KSP          ksp         
      Vec          bb
      PetscInt     ctx,n
      MatNullSpace nullspace
      PetscInt   xs,ys,zs,xm,ym,zm
      PetscErrorCode  ierr
      PetscScalar,pointer::rhs_array(:,:,:)
      PetscScalar,pointer::NullVecs_array(:)
      REAL :: dddwo,dddpo
      REAL :: ZP(IE,JE,KE),ZW(IE,JE,KEP)
      REAL :: XU(IE,JE),XP(IE,JE),YV(IE,JE),YP(IE,JE)  
                 
! pre process
      DO 2043 K=1,KEP
      DO 2043 J=1,JE
      DO 2043 I=1,IE
         ZW(I,J,K)=1.0/DZPNH(I,J,K)
2043  CONTINUE

      DO 2042 K=1,KE
      DO 2042 J=1,JE
      DO 2042 I=1,IE      	
         dddpo=(1.0/ZW(I,J,K) + 1.0/ZW(I,J,K+1))/2.0  
         ZP(I,J,K)=1.0/dddpo               	
2042  CONTINUE

      DO 2044 J=1,JE
      DO 2044 I=1,IE
        XU(I,J)=1.0/DLXU(I,J)
        XP(I,J)=1.0/DLXP(I,J)
        YV(I,J)=1.0/DLYV(I,J)
        YP(I,J)=1.0/DLYP(I,J)
2044  CONTINUE      

      call KSPGetDM(ksp,da,ierr)
      call DMDAGetCorners(da,xs,ys,zs,xm,ym,zm,ierr) 
      call DMDAVecGetArrayF90(da,bb,rhs_array,ierr)
      
      do 20 k=1,ke
      do 20 j=2,je1
      do 20 i=2,ie1
      	
#ifdef CG2D4NONHY           
        X_W=XU(I-1,J)*XP(I,J)
        X_C=(XU(I-1,J)+XU(I,J))*XP(I,J)
        X_E=XU(I,J)*XP(I,J)
        Y_N=YV(I,J-1)*YP(I,J)
        Y_C=(YV(I,J-1)+YV(I,J))*YP(I,J)
        Y_S=YV(I,J)*YP(I,J)
        Z_U=ZW(I,J,K)*ZP(I,J,K)
        Z_C=(ZW(I,J,K)+ZW(I,J,K+1))*ZP(I,J,K)
        Z_D=ZW(I,J,K+1)*ZP(I,J,K)      
      
        ! Centre 
        IF (WETO(I,J,K) .eq. 0) THEN  
         DIVG(I,J,K) = DIVG(I,J,K)+(X_C+Y_C+Z_C)*PNH_bot(I,J)   
        ELSE 
         ! west/east boundary
         IF     ( WETO(I-1,J,K).EQ.0 .AND. WETO(I+1,J,K).EQ.1 ) THEN
             X_C=XU(I,J)*XP(I,J)
          ELSEIF ( WETO(I-1,J,K).EQ.1 .AND. WETO(I+1,J,K).EQ.0 ) THEN
             X_C=XU(I-1,J)*XP(I,J)
          ELSEIF ( WETO(I-1,J,K).EQ.0 .AND. WETO(I+1,J,K).EQ.0 ) THEN
             X_C=0.0
          ENDIF
         ! north/south boundary
          IF     ( WETO(I,J-1,K).EQ.0 .AND. WETO(I,J+1,K).EQ.1 ) THEN
             Y_C=YV(I,J)*YP(I,J)
          ELSEIF ( WETO(I,J-1,K).EQ.1 .AND. WETO(I,J+1,K).EQ.0 ) THEN
             Y_C=YV(I,J-1)*YP(I,J)
          ELSEIF ( WETO(I,J-1,K).EQ.0 .AND. WETO(I,J+1,K).EQ.0 ) THEN
             Y_C=0.0
          ENDIF
         ! surface/bottom boundary
          IF (K .EQ. 1) THEN
             Z_U=0.0
          ELSEIF (K .EQ. KE) THEN
             Z_C=ZW(I,J,K)*ZP(I,J,K)
          ELSEIF (WETO(I,J,K+1).EQ.0) THEN
             Z_C=ZW(I,J,K)*ZP(I,J,K)
          ENDIF  
      
         ! Correction 
         ! Up     	     	    	     	      		     		
         IF (K .eq. 1 ) THEN
           DIVG(i,j,1) = DIVG(i,j,1) - Z_U*PNH_bot(I,J)
         elseif (WETO(I,J,K-1) .eq. 0) THEN
           DIVG(i,j,k) = DIVG(i,j,k) - Z_U*PNH_bot(I,J)
         ENDIF
         ! North      	
         IF (WETO(I,J-1,K) .eq. 0) THEN
           DIVG(i,j,k) = DIVG(i,j,k) - Y_N*PNH_bot(I,J-1)
         ENDIF
         ! West      	
         IF (WETO(I-1,J,K) .eq. 0) THEN
          DIVG(i,j,k) = DIVG(i,j,k) - X_W*PNH_bot(I-1,J)
         ENDIF

         ! East    	
         IF (WETO(I+1,J,K) .eq. 0) THEN
           DIVG(i,j,k) = DIVG(i,j,k) - X_E*PNH_bot(I+1,J)          
         ENDIF
         ! South        	
         IF (WETO(I,J+1,K) .eq. 0) THEN
           DIVG(i,j,k) = DIVG(i,j,k) - Y_S*PNH_bot(I,J+1)          
         ENDIF 
         ! Down           
         IF (K .eq. KE) THEN
           DIVG(i,j,k) = DIVG(i,j,k) - Z_D*PNH_bot(I,J)  
         elseif (WETO(I,J,K+1) .eq. 0) THEN 
           DIVG(i,j,k) = DIVG(i,j,k) - Z_D*PNH_bot(I,J)          
         ENDIF        
                   
        ENDIF
#endif /*CG2D4NONHY*/ 
       
      rhs_array(xs+i-2,ys+j-2,zs+k-1)=DIVG(i,j,k)*weto(i,j,k)    
      
20    continue
    
      call DMDAVecRestoreArrayF90(da,bb,rhs_array,ierr)  
         
      call VecAssemblyBegin(bb,ierr)
      call VecAssemblyEnd(bb,ierr)    
      
       if (Neumann)  then
         call VecCreateSeq(PETSC_COMM_SELF,ione,NullVecs,ierr)
         call VecGetArrayF90(NullVecs,NullVecs_array,ierr)
         NullVecs_array(ione)=zero_scalar
         call VecRestoreArrayF90(NullVecs,NullVecs_array,ierr)

         call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,izero,NullVecs,nullspace,ierr)
         call MatNullSpaceRemove(nullspace,bb,NullVecs,ierr)
         call MatNullSpaceDestroy(nullspace,ierr)
       endif
          
     return    
    END SUBROUTINE ComputeRHS



    SUBROUTINE compute_matrix_2d(ksp_2d,JJ,jac_2d,iopt,ierr)
      use petscmat
      DM   da_2d
      KSP  ksp_2d
      Mat  JJ,jac_2d
      PetscInt   n_petsc,xs,ys,xm,ym,mx,my
      PetscScalar  v(5)
      MatStencil   row(4),col(4,5)
      MatNullSpace nullspace
      PetscScalar,pointer::NullVecs_array(:)      
      REAL XU(IE,JE),XP(IE,JE),YV(IE,JE),YP(IE,JE)
                 
! pre process
      DO 2044 J=1,JE
      DO 2044 I=1,IE
        XU(I,J)=1.0/DLXU(I,J)
        XP(I,J)=1.0/DLXP(I,J)
        YV(I,J)=1.0/DLYV(I,J)
        YP(I,J)=1.0/DLYP(I,J)
2044  CONTINUE

! compute matrix
      call KSPGetDM(ksp_2d,da_2d,ierr)
      call DMDAGetInfo(da_2d,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,    &
      &                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                &
      &                 PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)      
      ! print*,'mx,my ',mx,my
      ! the corner indices : xs,ys,zs
      ! widths in the corresponding directions : xm,ym,zm       
      call DMDAGetCorners(da_2d,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
      ! print*,'xs,ys,xm,ym ',xs,ys,xm,ym
      DO 2046 I=2,IE1
      DO 2046 J=2,JE1
        X_W=XU(I-1,J)*XP(I,J)
        X_C=(XU(I-1,J)+XU(I,J))*XP(I,J)
        X_E=XU(I,J)*XP(I,J)
        Y_N=YV(I,J-1)*YP(I,J)
        Y_C=(YV(I,J-1)+YV(I,J))*YP(I,J)
        Y_S=YV(I,J)*YP(I,J)

        row(MatStencil_i) = xs+i-2
        row(MatStencil_j) = ys+j-2
        
       IF ( xs+i-2 .eq. 0 .OR. ys+j-2 .eq. 0 .OR. xs+i-2 .eq. mx-1 .OR. ys+j-2 .eq. my-1 ) THEN
          ! write to coefficient matrix on each Boundaries
          ! apply to Neumann Boundnary
                
          if (Neumann) then
             n_petsc=0
             X_CW=0
             X_CE=0
             Y_CN=0
             Y_CS=0
           if  ( ys+j-2 .ne. 0) then ! Non-North
               n_petsc = n_petsc + 1
               v(n_petsc) = Y_N
               Y_CN = YV(I,J-1)*YP(I,J)
               col(MatStencil_i,n_petsc) = row(MatStencil_i)
               col(MatStencil_j,n_petsc) = row(MatStencil_j)-1      
           endif
           if  ( xs+i-2 .ne. 0) then ! Non-West
                n_petsc = n_petsc + 1
                v(n_petsc) = X_W
                X_CW = XU(I-1,J)*XP(I,J)
                col(MatStencil_i,n_petsc) = row(MatStencil_i)-1
                col(MatStencil_j,n_petsc) = row(MatStencil_j)	
           endif       	
           if  ( xs+i-2 .ne. mx-1) then ! Non-East
               n_petsc = n_petsc + 1
               v(n_petsc) = X_E
               X_CE = XU(I,J)*XP(I,J)
               col(MatStencil_i,n_petsc) = row(MatStencil_i)+1
               col(MatStencil_j,n_petsc) = row(MatStencil_j)
           endif       
           if  (ys+j-2 .ne. my-1) then ! Non-South
               n_petsc = n_petsc + 1
               v(n_petsc) = Y_S
               Y_CS = YV(I,J)*YP(I,J)
              col(MatStencil_i,n_petsc) = row(MatStencil_i)
              col(MatStencil_j,n_petsc) = row(MatStencil_j)+1
           endif 
        
            n_petsc = n_petsc + 1
            v(n_petsc) = -(X_CW + X_CE + Y_CN + Y_CS)
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)
            call MatSetValuesStencil(jac_2d,1,row,n_petsc,col,v,INSERT_VALUES,ierr)
          else
          ! apply to Dirichlet Boundnary
            v(1) =-(X_C+Y_C)
            call MatSetValuesStencil(jac_2d,1,row,1,row,v,INSERT_VALUES,ierr)
          endif 
           
        ELSE
        ! write to coefficient matrix
          n_petsc=0
          ! west
            n_petsc=n_petsc+1
            v(n_petsc)     = X_W
            col(MatStencil_i,n_petsc) = row(MatStencil_i)-1
            col(MatStencil_j,n_petsc) = row(MatStencil_j)
          ! north
            n_petsc=n_petsc+1
            v(n_petsc)     = Y_N
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)-1
          ! center
            n_petsc=n_petsc+1
            v(n_petsc)     =-(X_C+Y_C)
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)
          ! south
            n_petsc=n_petsc+1
            v(n_petsc)     = Y_S
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)+1
          ! east
            n_petsc=n_petsc+1
            v(n_petsc)     = X_E
            col(MatStencil_i,n_petsc) = row(MatStencil_i)+1
            col(MatStencil_j,n_petsc) = row(MatStencil_j)              
          call MatSetValuesStencil(jac_2d,1,row,n_petsc,col,v,INSERT_VALUES,ierr)
          
        ENDIF
2046  CONTINUE
      call MatAssemblyBegin(jac_2d,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(jac_2d,MAT_FINAL_ASSEMBLY,ierr)
      
       if (Neumann) then
         call VecCreateSeq(PETSC_COMM_SELF,ione,NullVecs,ierr)
         call VecGetArrayF90(NullVecs,NullVecs_array,ierr)
            NullVecs_array(ione)=zero_scalar
         call VecRestoreArrayF90(NullVecs,NullVecs_array,ierr)   
         call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,izero,NullVecs,nullspace,ierr)
         call MatSetNullSpace(JJ,nullspace,ierr)
         call MatNullSpaceDestroy(nullspace,ierr)  
       endif  
           
      RETURN     
    END SUBROUTINE compute_matrix_2d 




    SUBROUTINE compute_matrix(ksp,JJ,jac,iopt,ierr)
      DM   da
      KSP  ksp
      Mat  JJ,jac
      PetscInt   n_petsc,xs,ys,zs,xm,ym,zm,mx,my,mz
      PetscScalar  v(7)
      MatStencil   row(4),col(4,7)
      MatNullSpace nullspace
      PetscScalar,pointer::NullVecs_array(:)         
      REAL ZP(IE,JE,KE),ZW(IE,JE,KEP)
      REAL XU(IE,JE),XP(IE,JE),YV(IE,JE),YP(IE,JE)
  
         
! pre process
      DO 2043 K=1,KEP
      DO 2043 J=1,JE
      DO 2043 I=1,IE
        ZW(I,J,K)=1.0/DZPNH(I,J,K)
2043  CONTINUE

      DO 2042 K=1,KE
      DO 2042 J=1,JE
      DO 2042 I=1,IE      
         dddpo=(1.0/ZW(I,J,K) + 1.0/ZW(I,J,K+1))/2.0  
         ZP(I,J,K)=1.0/dddpo               
2042  CONTINUE


      DO 2044 J=1,JE
      DO 2044 I=1,IE
        XU(I,J)=1.0/DLXU(I,J)
        XP(I,J)=1.0/DLXP(I,J)
        YV(I,J)=1.0/DLYV(I,J)
        YP(I,J)=1.0/DLYP(I,J)
2044  CONTINUE

      call KSPGetDM(ksp,da,ierr)   
      call DMDAGetCorners(da,xs,ys,zs,xm,ym,zm,ierr)
!      print*,xs,ys,zs,xm,ym,zm
      call DMDAGetInfo(da,PETSC_NULL_INTEGER,mx,my,mz,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
      & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
      & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
      DO 2046 I=2,IE1
      DO 2046 J=2,JE1
      DO 2046 K=1,KE
        X_W=XU(I-1,J)*XP(I,J)
        X_C=(XU(I-1,J)+XU(I,J))*XP(I,J)
        X_E=XU(I,J)*XP(I,J)
        Y_N=YV(I,J-1)*YP(I,J)
        Y_C=(YV(I,J-1)+YV(I,J))*YP(I,J)
        Y_S=YV(I,J)*YP(I,J)
        Z_U=ZW(I,J,K)*ZP(I,J,K)
        Z_C=(ZW(I,J,K)+ZW(I,J,K+1))*ZP(I,J,K)
        Z_D=ZW(I,J,K+1)*ZP(I,J,K)

        row(MatStencil_i) = xs+i-2
        row(MatStencil_j) = ys+j-2
        row(MatStencil_k) = zs+k-1
        IF (weto(I,J,K) .eq. 0 ) THEN
! The land ponits do not need to discretized.        
          v(1) =-(X_C+Y_C+Z_C)
          call MatSetValuesStencil(jac,1,row,1,row,v,INSERT_VALUES,ierr)
        ELSE
! The homogeneous Neumann Boundary Conditions are accepted for the land boundaries
          ! west/east boundary
          IF     ( weto(I-1,J,K).EQ.0 .AND. weto(I+1,J,K).EQ.1 ) THEN
             X_W=0.0
             X_C=XU(I,J)*XP(I,J)
          ELSEIF ( weto(I-1,J,K).EQ.1 .AND. weto(I+1,J,K).EQ.0 ) THEN
             X_E=0.0
             X_C=XU(I-1,J)*XP(I,J)
          ELSEIF ( WETO(I-1,J,K).EQ.0 .AND. WETO(I+1,J,K).EQ.0 ) THEN
             X_W=0.0
             X_E=0.0
             X_C=0.0
          ENDIF
          
          ! north/south boundary          
          IF     ( WETO(I,J-1,K).EQ.0 .AND. WETO(I,J+1,K).EQ.1 ) THEN
             Y_N=0.0
             Y_C=YV(I,J)*YP(I,J)
          ELSEIF ( WETO(I,J-1,K).EQ.1 .AND. WETO(I,J+1,K).EQ.0 ) THEN
             Y_S=0.0
             Y_C=YV(I,J-1)*YP(I,J)
          ELSEIF ( WETO(I,J-1,K).EQ.0 .AND. WETO(I,J+1,K).EQ.0 ) THEN
             Y_N=0.0
             Y_S=0.0
             Y_C=0.0
          ENDIF

! north/south boundary for purely 2-D (ie_g*1*ke)
          !if (ys+j-2 .eq. 0 ) then ! north
          !    Y_N=0             
          !    Y_S=Y_C              
          !endif
                
          !if (ys+j-2 .eq. my-1 ) then ! south  
          !    Y_N=Y_C 
          !    Y_S=0      
          !endif
          
          ! surface/bottom boundary
          IF (K .EQ. 1) THEN
              Z_U=0.0
          ELSEIF (K .EQ. KE) THEN
             Z_D=0.0
             Z_C=ZW(I,J,K)*ZP(I,J,K)
          ELSEIF (weto(I,J,K+1).EQ.0) THEN
             Z_D=0.0
             Z_C=ZW(I,J,K)*ZP(I,J,K)
          ENDIF

! No treating for topo, dpnh/dx=0 or pnh=0 at OB (homogeneous Neumann or Dirichlet Conditions)
! Open Boundary Conditions for Solving Poission Problem
! It also refers to the boundary of region for sloving Poisson Problem          
          if (ICYCLI_X.EQ.0) then
            if (have_g_is .AND. West .AND. (I.EQ.2)) then
              X_W=0.0
!              X_C=XU(I,J)*XP(I,J)
            endif
            if (have_g_ie .AND. East .AND. (I.EQ.IE1)) then
              X_E=0.0
!              X_C=XU(I-1,J)*XP(I,J)
             endif
          endif
          
          if (ICYCLI_Y.EQ.0) then
            if (have_g_js .AND. North .AND. (J.EQ.2)) then
              Y_N=0.0
!              Y_C=YV(I,J)*YP(I,J)
            endif
            if (have_g_je .AND. South .AND. (J.EQ.JE1)) then
              Y_S=0.0
!              Y_C=YV(I,J-1)*YP(I,J)
            endif
          endif

   ! write to coefficient matrix
          n_petsc=0
          ! west
          IF (weto(I-1,J,K) .EQ. 1 &
#ifndef CYCLICWE          
             & .AND. xs+i-2 .NE. 0   &
#endif          
             & ) THEN
            n_petsc=n_petsc+1
            v(n_petsc)     = X_W
            col(MatStencil_i,n_petsc) = row(MatStencil_i)-1
            col(MatStencil_j,n_petsc) = row(MatStencil_j)
            col(MatStencil_k,n_petsc) = row(MatStencil_k)
          ENDIF
          ! north
          IF (weto(I,J-1,K) .EQ. 1 &
#ifndef CYCLICNS           
             & .AND. ys+j-2 .NE. 0 &
#endif          
             & ) THEN
            n_petsc=n_petsc+1
            v(n_petsc)     = Y_N
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)-1
            col(MatStencil_k,n_petsc) = row(MatStencil_k)
          ENDIF
          ! up
          IF (K .NE. 1) THEN
            n_petsc=n_petsc+1
            v(n_petsc)     = Z_U
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)
            col(MatStencil_k,n_petsc) = row(MatStencil_k)-1
          ENDIF
          ! center
            n_petsc=n_petsc+1
            v(n_petsc)     =-(X_C+Y_C+Z_C)
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)
            col(MatStencil_k,n_petsc) = row(MatStencil_k)
          ! down
          IF (K .NE. KE) THEN
          IF (weto(I,J,K+1) .EQ. 1) THEN
            n_petsc=n_petsc+1
            v(n_petsc)     = Z_D
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)
            col(MatStencil_k,n_petsc) = row(MatStencil_k)+1
          ENDIF
          ENDIF
          ! South
          IF (weto(I,J+1,K) .EQ. 1 &
#ifndef CYCLICNS            
             & .AND. ys+j-2 .NE. my-1 &
#endif          
             & ) THEN
            n_petsc=n_petsc+1
            v(n_petsc)     = Y_S
            col(MatStencil_i,n_petsc) = row(MatStencil_i)
            col(MatStencil_j,n_petsc) = row(MatStencil_j)+1
            col(MatStencil_k,n_petsc) = row(MatStencil_k)
          ENDIF
          ! east
          IF (weto(I+1,J,K) .EQ. 1 &
#ifndef CYCLICWE            
             & .AND. xs+i-2 .NE. mx-1 &
#endif           
             & ) THEN
            n_petsc=n_petsc+1
            v(n_petsc)     = X_E
            col(MatStencil_i,n_petsc) = row(MatStencil_i)+1
            col(MatStencil_j,n_petsc) = row(MatStencil_j)
            col(MatStencil_k,n_petsc) = row(MatStencil_k)
          ENDIF
          call MatSetValuesStencil(jac,1,row,n_petsc,col,v,INSERT_VALUES,ierr)
        ENDIF
2046  CONTINUE
      call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
      
       if (Neumann) then
         call VecCreateSeq(PETSC_COMM_SELF,ione,NullVecs,ierr)
         call VecGetArrayF90(NullVecs,NullVecs_array,ierr)
            NullVecs_array(ione)=zero_scalar
         call VecRestoreArrayF90(NullVecs,NullVecs_array,ierr)   
         call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,izero,NullVecs,nullspace,ierr)
         call MatSetNullSpace(JJ,nullspace,ierr)
         call MatNullSpaceDestroy(nullspace,ierr)  
       endif        
      
      RETURN
    END SUBROUTINE compute_matrix

  

    SUBROUTINE MyKSPConverged(ksp,n,rnorm,flag,dummy,ierr)
    IMPLICIT NONE
!  Input Parameters:
!    ksp   - iterative context
!    n     - iteration number
!    rnorm - 2-norm (preconditioned) residual value (may be estimated)
!    dummy - optional user-defined monitor context (unused here)
    KSP ksp
    PetscErrorCode ierr
    PetscInt n,dummy
    PetscReal rnorm
    KSPConvergedReason flag

    if (rnorm .le. .05) then
         flag = 1
    else
         flag = 0
    endif
    write(IO_STDOUT,*) 'MyKSPConvergedTest::'
    write(IO_STDOUT,*) '  Norm of error:', rnorm
    write(IO_STDOUT,*) 'ConvergedReason:',flag
    write(IO_STDOUT,*) 'IterationNumber:', n   
       
    END SUBROUTINE MyKSPConverged
#endif
END MODULE MO_PETSC

