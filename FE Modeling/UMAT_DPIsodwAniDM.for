C     This UMAT subroutine performs the incremental integration of elastoplastic damage model with elastic stiffness dependent on the stress tensor (AniDM)
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
!     Implicit real*8(a-h,o-z)
!     This means that, by default, any variables with
!     first letter between a-h or o-z are double precision.
!     The rest are integers.
C
      CHARACTER*80 CMNAME
      REAL*8 NORM_SEFF_TR,LAMBDA_P,NORM_EPSPLUS,FV1(3),EIGVALR(3),
     2 EIGVALI(3),EIGVEC(3,3),CKAPPA_D,LAMBDA_KM1,LAMBDA_K,LAMBDA_KP1,
     3 LAMBDA,NORM_SEFF,I1,NU0,NORM_RESID,LAMBDA2,I2I2(3,3,3,3),
     4 SIGEFF_OLD(3,3)
      !DDSDDT,DRPLED,
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      !JSTEP(4)
C	  
      DIMENSION EPS_OLD(3,3),DEPS(3,3),EPS(3,3),
     1 SIGEFF_TR(3,3),SEFF_TR(3,3),SIGEFF(3,3),SEFF(3,3),
     2 C1(3,3,3,3),DSIG(3,3),DEPS_E(3,3),DEPS_P(3,3),AHAT(3,3),
     3 SIG_OLD(3,3),DSIG_TR(3,3),ZERO_T(3,3),C_EFF0(NTENS,NTENS),
     4 SIGEFF_OLD_V(NTENS),SIG(3,3),SIGEFF_V(NTENS),C3333(3,3,3,3),
     5 C_EFF(NTENS,NTENS),EYE2(3,3),DSTRESS(NTENS),STRESS_OLD(NTENS),
     6 PPLUS(3,3,3,3),EPSPLUS(3,3),SIG_TR(3,3),IV1(3),EPS_E_OLD(3,3),
     7 C_SD(3,3,3,3),CNN(NTENS,NTENS),EPS_E_OLD_V(NTENS),EPS_E(3,3),
     8 C(3,3,3,3),EPS_E_V(NTENS),CR(3,3,3,3),PC1(3,3,3,3),PC1P(3,3,3,3),
     9 EPS_P(3,3),EPS_P_OLD(3,3),SEFFSEFF(3,3,3,3),DF_DSIG(3,3),
     1 DF2_DSIG2(3,3,3,3),DF2_DSIG266(6,6),S_SD66(6,6),
     2 XI(3,3,3,3),XIDF_DSIG(3,3),DEL_SIGEFF(3,3),RESID(3,3),R(3),
     3 R_tr(3),P_dev(3,3,3,3),Rl2d(3,3),Qij(3,3),Pij(3,3,3,3),
     4 EYE4(3,3,3,3),C_SD66(6,6),S_SD(3,3,3,3),DEL_EPS_P(3,3),
     5 STRAN_E(NTENS),ALPHA(3),ALPHA_TR(3),P3333(3,3,3,3),Q_STAR(3,3),
     6 Q_REG(3,3),Q_K(3,3),XI_INV66(6,6),XI66(6,6)
  
      COMMON/PARAMETERS/ELAMBDA0,EMU0,CKAPPA,C0D,C1D,C2D,B,R0,R1,R2,EK0
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
C
      PARAMETER (PI=3.1415926526D0,nwrite=1)
C
C     MATERIAL PARAMETERS
C
      !WRITE(7,*) 'NTENS=',NTENS
      E0 = PROPS(1) !E0 is E1 in power law formulation of stress dependency
      NU0 = PROPS(2)
      CKAPPA = PROPS(3)
      C0D = PROPS(4) 
      C1D = PROPS(5)
      C2D = PROPS(6)
      bb   = PROPS(7)
      B   = PROPS(8)
      
      ELAMBDA0 = E0*NU0/((1.+NU0)*(1.-2.*NU0))
      EMU0 = E0/(2.*(1.+NU0))
      EK0 = ELAMBDA0+2./3.*EMU0

C     RPL=1
      DO I=1,3
        DO J=1,3
          ZERO_T(I,J)=0.0D0
        END DO
      END DO
      CALL IDENTITYMAT2(3,EYE2)
C
C====================================================================== 
C     DEFAULT NTENS = 6 FOR 3D CASE, NTENS = 4 FOR FOR PLANE STRAIN CASES
      ! *_V: vectorial form
      DO I = 1, NTENS   
        SIGEFF_OLD_V(I) = STATEV(I)       !Recover effective stress
        !STRESS(I) = STATEV(I+3+NTENS)         !Recover true stress 
        EPS_E_OLD_V(I)=STATEV(NTENS+I)
      END DO
      D = STATEV(1+2*NTENS)
      DO I = 1,3
          ALPHA(I)=STATEV(I+2*NTENS+5)   !Store the effective stress
      END DO
      
      IF (TIME(2).EQ.ZERO) THEN
          ALPHA(:)=ONE ! ASSUME ALPHA=1 AT BEGINING
          D=ZERO ! NO DAMAGE AT BEGINING
      END IF
      
      ! Avoid any floating error in stress tensor
      DO I=1,NTENS
          IF (ABS(SIGEFF_OLD_V(I)).LT.1E-8) THEN
              SIGEFF_OLD_V(I) =ZERO
          ENDIF
      ENDDO
      
      !CONVERT VECTOR TO 2ND-RANK TENSOR
      CALL MAT1_MAT2(NTENS,SIGEFF_OLD_V,SIGEFF_OLD,ONE)
      CALL MAT1_MAT2(NTENS,EPS_E_OLD_V,EPS_E_OLD,ONE)
      CALL MAT1_MAT2(NTENS,DSTRAN,DEPS,HALF)
C     write(7,*) 'DEPS=',DEPS
      CALL MAT1_MAT2(NTENS,STRAN,EPS_OLD,HALF)
      CALL MAT1_MAT2(NTENS,STRESS,SIG_OLD,ONE)
C
C     UPDATE DAMAGE VARIABLE
      EPS=EPS_OLD+DEPS
      D_OLD=D
      CALL D_FCN(EPS,D,FD)
      CKAPPA_D=(1-D)*CKAPPA
C      
C --- CALCULATION OF EFFECTIVE FOURTH-ORDER ELASTIC STIFFNESS
C     TENSOR AND 2ND-ORDER STRESS TENSOR WITH DAMAGE EXCLUDED
      CALL CAL_C(NTENS,ZERO,ZERO,C1,C_EFF0)
C     CALCULATE DAMAGE EFFECTIVE STRESS (TRIAL)
      EPS_E=EPS_E_OLD+DEPS
      
C_____(Activate iteration to update alpha in each step)_____
      DEL_ALPHA=1
      TOL_ALPHA=1E-5
      DO WHILE (ABS(DEL_ALPHA)>TOL_ALPHA)
C___________________________________________________________
          ALPHA_TR=ALPHA
          CALL GEIGEN(3,SIGEFF_OLD,EIGVALR,EIGVALI,EIGVEC,FV1,IV1)
          !EIGVEC=EYE2
          Q_REG=ZERO
          Q_STAR=ZERO
          DO K=1,3
              CALL Ai_Bj(3,EIGVEC(:,K),EIGVEC(:,K),Q_K)
              Q_REG=Q_REG+Q_K
              Q_STAR=Q_STAR+ALPHA_TR(K)**(1./4.)*Q_K
          ENDDO
          
          DO I = 1, 3   
              DO J=1, 3
                  DO K=1,3
                      DO L=1,3
                          P3333(I,J,K,L)=ZERO
                          DO M=1,3
                              DO N=1,3
                          P3333(I,J,K,L)=P3333(I,J,K,L)+
     2                    Q_STAR(I,M)*Q_STAR(J,N)*Q_REG(K,M)*Q_REG(L,N)            
                              ENDDO
                          ENDDO
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
          CALL Aijkl_Bklmn(3,P3333,C1,PC1)
          CALL Aijkl_Bklmn(3,PC1,P3333,C_SD)
          CALL Aijkl_Bkl(3,C_SD,EPS_E,SIGEFF_TR)
          CALL TRACE(3,SIGEFF_TR,SIGEFF_TRii)
          p_tr=SIGEFF_TRii/3
          SEFF_TR=SIGEFF_TR-p_tr*EYE2

          CALL Aij_Bij(3,SEFF_TR,SEFF_TR,SEFF_TR2)
          NORM_SEFF_TR=SQRT(SEFF_TR2)
          q=SQRT(3./2.)*NORM_SEFF_TR
C         EVALUATE YIELD FUNCTION
          !IF (p_tr.GE.0) THEN !yield in tension
          !    FP=q-CKAPPA_D
          !ELSE
          !    FP=q+B*p_tr-CKAPPA_D  !yield in compression
          !ENDIF
          FP=q+B*p_tr-CKAPPA_D 
          EPS_P_OLD=EPS_OLD-EPS_E_OLD
C         UPDATE DAMAGE EFFECTIVE STRESS
          IF (FP.LT.0) THEN !NO YIELDING
              SIGEFF=SIGEFF_TR
              EPS_P=EPS_P_OLD
              LAMBDA=0
              p=p_tr
  
          ELSE !YIELDING
              ! Computation constants
              LAMBDA=0
              CALL IDENTITYMAT4(3,EYE4)
              CALL Aij_Bkl(3,EYE2,EYE2,I2I2)
              P_dev=EYE4-1/3*I2I2
              TOL_FP=1E-8
              TOL_RESID=1E-8
              CALL MAT4_MAT2(6,C_SD,C_SD66,1)
              CALL INVERSE(C_SD66,6,6,S_SD66)
              CALL MAT2_MAT4(6,S_SD66,S_SD,2)
              
              ! Initialize the iteration
              EPS_P=EPS_P_OLD
              p=p_tr
              SEFF=SEFF_TR
              SIGEFF=SEFF+p*EYE2
              
              CALL Aij_Bkl(3,SEFF,SEFF,SEFFSEFF)
              CALL Aij_Bij(3,SEFF_TR,SEFF_TR,SEFF_TR2)
              NORM_SEFF_TR=SQRT(SEFF_TR2)
              AHAT=SEFF/NORM_SEFF_TR
              DF_DSIG=SQRT(3./2.)*AHAT+B/3*EYE2
              DF2_DSIG2=SQRT(3./2.)*(P_dev/NORM_SEFF_TR-
     2        SEFFSEFF/NORM_SEFF_TR**3)
              RESID=-EPS_P+EPS_P_OLD+LAMBDA*DF_DSIG
              CALL Aij_Bij(3,RESID,RESID,RESID2)
              NORM_RESID=SQRT(RESID2)
              
              ! Debugging Break
              !dbgVar=1
              !Do While(dbgVar .ne. 999)
              !    dbgVar = 1
              !End Do
              
              !SIMO & HUGHES 1997
              DO WHILE ((ABS(FP).GT.TOL_FP).OR.(LAMBDA.LT.0).OR.
     2        (NORM_RESID.GT.TOL_RESID))
                  
                  ! STEP 1
                  CALL MAT4_MAT2(6,DF2_DSIG2,DF2_DSIG266,2)
                  XI_INV66=S_SD66+LAMBDA*DF2_DSIG266
                  CALL INVERSE(XI_INV66,6,6,XI66)
                  CALL MAT2_MAT4(6,XI66,XI,1)
                  CALL Aijkl_Bkl(3,XI,DF_DSIG,XIDF_DSIG)
                  CALL Aij_Bij(3,RESID,XIDF_DSIG,RXD)
                  CALL Aij_Bij(3,DF_DSIG,XIDF_DSIG,DXD)
                  LAMBDA2=(FP-RXD)/DXD
                  
                  ! STEP 2
                  Rl2d=-RESID-LAMBDA2*DF_DSIG
                  CALL Aijkl_Bkl(3,XI,Rl2d,DEL_SIGEFF)
                  CALL Aijkl_Bkl(3,-S_SD,DEL_SIGEFF,DEL_EPS_P)
                  
                  ! STEP 3
                  EPS_P=EPS_P+DEL_EPS_P
                  LAMBDA=LAMBDA+LAMBDA2
                  SIGEFF=SIGEFF+DEL_SIGEFF
                  
                  ! Evaluate termination criteria
                  CALL TRACE(3,SIGEFF,SIGEFFii)
                  p=SIGEFFii/3
                  SEFF=SIGEFF-p*EYE2
                  CALL Aij_Bkl(3,SEFF,SEFF,SEFFSEFF)
                  CALL Aij_Bij(3,SEFF,SEFF,SEFF2)
                  NORM_SEFF=SQRT(SEFF2)
                  AHAT=SEFF/NORM_SEFF
                  DF_DSIG=SQRT(3./2.)*AHAT+B/3*EYE2
                  DF2_DSIG2=SQRT(3./2.)*(P_dev/NORM_SEFF-
     2            SEFFSEFF/NORM_SEFF**3)
                  FP=SQRT(3./2.)*NORM_SEFF+B*p-CKAPPA_D
                  RESID=-EPS_P+EPS_P_OLD+LAMBDA*DF_DSIG
                  CALL Aij_Bij(3,RESID,RESID,RESID2)
                  NORM_RESID=SQRT(RESID2)
              ENDDO
      
          !ENDIF
          ENDIF
      
C     UPDATE STRESS-DEPENDENT STIFFNESS     
          EPS_E=EPS-EPS_P
          CALL GEIGEN(3,SIGEFF,EIGVALR,EIGVALI,EIGVEC,FV1,IV1)
          DO K=1,3
              SIGEFF_K=EIGVALR(K)
              IF (SIGEFF_K>-0.01) THEN
                  SIGEFF_K=-0.01
              ENDIF
              ALPHA_K=(-SIGEFF_K)**bb
              !ALPHA_K=1+(-SIGEFF_K)*bb
              ALPHA(K)=ALPHA_K
          END DO
C_____(Activate iteration to update alpha in each step)_____
          DEL_ALPHA=NORM2(ALPHA-ALPHA_TR)
          
          !THE FOLLOWING LINE IS TO SKIP ITERATION (COMMENT IT TO ENABLE ITERATION)
          DEL_ALPHA=ZERO 
          
      ENDDO
C___________________________________________________________
      
      C=C_SD*(1-D)
      CALL Aijkl_Bkl(3,C,EPS_E,SIG)
      CALL MAT4_MAT2(NTENS,C,CNN,1) 
      
C
C --- UPDATING STATE VARIABLES
C     
      CALL MAT2_MAT1(NTENS,SIG,STRESS,ONE)
      CALL MAT2_MAT1(NTENS,EPS_E,EPS_E_V,ONE)
      CALL MAT2_MAT1(NTENS,SIGEFF,SIGEFF_V,ONE)
      CALL MAT2_MAT1(NTENS,EPS_E,STRAN_E,TWO)
      !CALL Aij_Bj(NTENS,CNN,STRAN_E,STRESS)
C     
C  STORE STATE VARIABLES
C
      DO I = 1,NTENS
        STATEV(I) = SIGEFF_V(I)   !Store the effective stress
        STATEV(I+NTENS) = EPS_E_V(I)   !
      END DO
      STATEV(1+2*NTENS) = D
      STATEV(2+2*NTENS) = LAMBDA

      DO I = 1,3
          STATEV(I+2*NTENS+5) = ALPHA(I)   ! ALPHA(1)=STATEV(14) to avoid STATEV(12) (Abaqus Bug)
      END DO
      STATEV(2*NTENS+9) = (ALPHA(1)+ALPHA(2)+ALPHA(3))/3
      
C --- UPDATING THE JACOBIAN MATRIX FOR ACCELERATING CONVERGENCE
      !Use the current stiffness as the Jacobian
      DDSDDE=CNN/(1-bb)

      IF (NOEL.EQ.1) THEN !ASCE Wall
          IF (NPT.EQ.3) THEN
              WRITE(7,*) 'DDSDDE='
              DO i=1,NTENS
                  WRITE(7,"(100g15.8)") (DDSDDE(i,j), j=1,NTENS)
              ENDDO
              WRITE(7,*) 'C_EFF='
              DO i=1,NTENS
                  WRITE(7,"(100g15.8)") (CNN(i,j), j=1,NTENS)
              ENDDO
              !WRITE(7,*) 'CD_ACT='
              !DO i=1,NTENS
              !    WRITE(7,"(100g15.8)") (CD_ACT_NN(i,j), j=1,NTENS)
              !ENDDO
              WRITE(7,*) 'EPS='
              DO i=1,3
                  WRITE(7,"(100g15.8)") (EPS(i,j), j=1,3)
              ENDDO
              WRITE(7,*) 'p='
              WRITE(7,"(100g15.8)") p
              WRITE(7,*) 'D='
              WRITE(7,"(100g15.8)") D
              WRITE(7,*) 'alpha='
              DO i=1,3
                  WRITE(7,"(100g15.8)") (ALPHA(i))
              ENDDO
              WRITE(7,*) 'FD='
              WRITE(7,"(100g15.8)") FD
              WRITE(7,*) 'FP='
              WRITE(7,"(100g15.8)") FP  
              WRITE(7,*) 'SIGEFF_OLD='
              DO i=1,3
                  WRITE(7,"(100g15.8)") (SIGEFF_OLD(i,j), j=1,3)
              ENDDO
              WRITE(7,*) 'EIGVEC='
              DO i=1,3
                  WRITE(7,"(100g15.8)") (EIGVEC(i,j), j=1,3)
              ENDDO
          ENDIF
      ENDIF

      RETURN
      END
C
      SUBROUTINE MAT1_MAT2(NTENS,VECTOR,TENSOR,FACT)
C ====================================================================
C 
C                MAT1_MAT2 : VECTOR TO TENSOR  
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
        TENSOR=0.0D0
C
        TENSOR(1 , 1) = VECTOR( 1 )
        TENSOR(2 , 2) = VECTOR( 2 )
        TENSOR(3 , 3) = VECTOR( 3 )
        TENSOR(1 , 2) = VECTOR( 4 )*FACT
        TENSOR(2 , 1) = TENSOR(1 , 2)
        IF (NTENS>4) THEN
            TENSOR(1 , 3) = VECTOR( 5 )*FACT
            TENSOR(3 , 1) = TENSOR(1 , 3)
            TENSOR(2 , 3) = VECTOR( 6 )*FACT
            TENSOR(3 , 2) = TENSOR(2 , 3)
        ENDIF
C
      RETURN
      END
C
      SUBROUTINE MAT4_MAT2(NTENS,TENSOR,DMATRIX,ICOE)
C
C ====================================================================
C                        MAT4_MAT2                                                                  I
C I        THIS PROGRAM TRANSFORMS THE FOURTH ORDER COMPLIANCE       I
C I        TENSOR TO A SECOND ORDER MATRIX                           I
C I                                                                  I
C ====================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TENSOR(3,3,3,3),DMATRIX(NTENS,NTENS)
C
      DATA ZERO,TWO /0.0D0,2.0D0/
C
C     D2 = THE SECOND ORDER STIFFNESS MATRIX
C
      DMATRIX=0.0D0

      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C
      DMATRIX(1,1)=TENSOR(1,1,1,1)
      DMATRIX(1,2)=TENSOR(1,1,2,2)
      DMATRIX(1,3)=TENSOR(1,1,3,3)
      DMATRIX(1,4)=TENSOR(1,1,1,2)*COE1
      DMATRIX(2,1)=TENSOR(2,2,1,1)
      DMATRIX(2,2)=TENSOR(2,2,2,2)
      DMATRIX(2,3)=TENSOR(2,2,3,3)
      DMATRIX(2,4)=TENSOR(2,2,1,2)*COE1
      DMATRIX(3,1)=TENSOR(3,3,1,1)
      DMATRIX(3,2)=TENSOR(3,3,2,2)
      DMATRIX(3,3)=TENSOR(3,3,3,3)
      DMATRIX(3,4)=TENSOR(3,3,1,2)*COE1
      DMATRIX(4,1)=TENSOR(1,2,1,1)*COE1
      DMATRIX(4,2)=TENSOR(1,2,2,2)*COE1
      DMATRIX(4,3)=TENSOR(1,2,3,3)*COE1
      DMATRIX(4,4)=TENSOR(1,2,1,2)*COE2
C
      IF (NTENS>4) THEN
          DMATRIX(1,5)=TENSOR(1,1,2,3)*COE1
          DMATRIX(1,6)=TENSOR(1,1,1,3)*COE1
          DMATRIX(2,5)=TENSOR(2,2,2,3)*COE1
          DMATRIX(2,6)=TENSOR(2,2,1,3)*COE1
          DMATRIX(3,5)=TENSOR(3,3,2,3)*COE1
          DMATRIX(3,6)=TENSOR(3,3,1,3)*COE1
          DMATRIX(4,5)=TENSOR(1,2,2,3)*COE2
          DMATRIX(4,6)=TENSOR(1,2,1,3)*COE2
C
          DMATRIX(5,1)=TENSOR(2,3,1,1)*COE1
          DMATRIX(5,2)=TENSOR(2,3,2,2)*COE1
          DMATRIX(5,3)=TENSOR(2,3,3,3)*COE1
          DMATRIX(5,4)=TENSOR(2,3,1,2)*COE2
          DMATRIX(5,5)=TENSOR(2,3,2,3)*COE2
          DMATRIX(5,6)=TENSOR(2,3,1,3)*COE2
C
          DMATRIX(6,1)=TENSOR(1,3,1,1)*COE1
          DMATRIX(6,2)=TENSOR(1,3,2,2)*COE1
          DMATRIX(6,3)=TENSOR(1,3,3,3)*COE1
          DMATRIX(6,4)=TENSOR(1,3,1,2)*COE2
          DMATRIX(6,5)=TENSOR(1,3,2,3)*COE2
          DMATRIX(6,6)=TENSOR(1,3,1,3)*COE2
      ENDIF
C
C
      RETURN
      END
C
      SUBROUTINE MAT2_MAT4(NTENS,DMATRIX,TENSOR,ICOE)
C======================================================================================
C
C                             MAT2_MAT4
C
C======================================================================================
      INCLUDE 'ABA_PARAM.INC'
C  
      DIMENSION DMATRIX(6,6),TENSOR(3,3,3,3)
C
C     INITALIZATION
C
       DO I=1,3
         DO J=1,3
           DO K=1,3
             DO L=1,3 
               TENSOR(I,J,K,L)=0.0D0
             END DO
           END DO 
         END DO
       END DO
C
      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C      
      TENSOR(1,1,1,1) = DMATRIX(1,1)
      TENSOR(1,1,2,2) = DMATRIX(1,2)
      TENSOR(1,1,3,3) = DMATRIX(1,3)
      TENSOR(1,1,1,2) = DMATRIX(1,4)/COE1
      TENSOR(1,1,2,1) = DMATRIX(1,4)/COE1
      TENSOR(1,1,2,3) = DMATRIX(1,5)/COE1
      TENSOR(1,1,3,2) = DMATRIX(1,5)/COE1
      TENSOR(1,1,1,3) = DMATRIX(1,6)/COE1
      TENSOR(1,1,3,1) = DMATRIX(1,6)/COE1
C
      TENSOR(2,2,1,1) = DMATRIX(2,1)
      TENSOR(2,2,2,2) = DMATRIX(2,2)
      TENSOR(2,2,3,3) = DMATRIX(2,3)
      TENSOR(2,2,1,2) = DMATRIX(2,4)/COE1
      TENSOR(2,2,2,1) = DMATRIX(2,4)/COE1
      TENSOR(2,2,2,3) = DMATRIX(2,5)/COE1
      TENSOR(2,2,3,2) = DMATRIX(2,5)/COE1
      TENSOR(2,2,1,3) = DMATRIX(2,6)/COE1
      TENSOR(2,2,3,1) = DMATRIX(2,6)/COE1
C
      TENSOR(3,3,1,1) = DMATRIX(3,1)
      TENSOR(3,3,2,2) = DMATRIX(3,2)
      TENSOR(3,3,3,3) = DMATRIX(3,3)
      TENSOR(3,3,1,2) = DMATRIX(3,4)/COE1
      TENSOR(3,3,2,1) = DMATRIX(3,4)/COE1
      TENSOR(3,3,2,3) = DMATRIX(3,5)/COE1
      TENSOR(3,3,3,2) = DMATRIX(3,5)/COE1
      TENSOR(3,3,1,3) = DMATRIX(3,6)/COE1
      TENSOR(3,3,3,1) = DMATRIX(3,6)/COE1
C
      TENSOR(1,2,1,1) = DMATRIX(4,1)/COE1
      TENSOR(1,2,2,2) = DMATRIX(4,2)/COE1
      TENSOR(1,2,3,3) = DMATRIX(4,3)/COE1
      TENSOR(1,2,1,2) = DMATRIX(4,4)/COE2
      TENSOR(1,2,2,1) = DMATRIX(4,4)/COE2
      TENSOR(1,2,2,3) = DMATRIX(4,5)/COE2
      TENSOR(1,2,3,2) = DMATRIX(4,5)/COE2
      TENSOR(1,2,1,3) = DMATRIX(4,6)/COE2
      TENSOR(1,2,3,1) = DMATRIX(4,6)/COE2
C
      TENSOR(2,3,1,1) = DMATRIX(5,1)/COE1
      TENSOR(2,3,2,2) = DMATRIX(5,2)/COE1
      TENSOR(2,3,3,3) = DMATRIX(5,3)/COE1
      TENSOR(2,3,1,2) = DMATRIX(5,4)/COE2
      TENSOR(2,3,2,1) = DMATRIX(5,4)/COE2
      TENSOR(2,3,2,3) = DMATRIX(5,5)/COE2
      TENSOR(2,3,3,2) = DMATRIX(5,5)/COE2
      TENSOR(2,3,1,3) = DMATRIX(5,6)/COE2
      TENSOR(2,3,3,1) = DMATRIX(5,6)/COE2
C
      TENSOR(1,3,1,1) = DMATRIX(6,1)/COE1
      TENSOR(1,3,2,2) = DMATRIX(6,2)/COE1
      TENSOR(1,3,3,3) = DMATRIX(6,3)/COE1
      TENSOR(1,3,1,2) = DMATRIX(6,4)/COE2
      TENSOR(1,3,2,1) = DMATRIX(6,4)/COE2
      TENSOR(1,3,2,3) = DMATRIX(6,5)/COE2
      TENSOR(1,3,3,2) = DMATRIX(6,5)/COE2
      TENSOR(1,3,1,3) = DMATRIX(6,6)/COE2
      TENSOR(1,3,3,1) = DMATRIX(6,6)/COE2
C      
      TENSOR(2,1,1,1) = DMATRIX(4,1)/COE1
      TENSOR(2,1,2,2) = DMATRIX(4,2)/COE1
      TENSOR(2,1,3,3) = DMATRIX(4,3)/COE1
      TENSOR(2,1,1,2) = DMATRIX(4,4)/COE2
      TENSOR(2,1,2,1) = DMATRIX(4,4)/COE2
      TENSOR(2,1,2,3) = DMATRIX(4,5)/COE2
      TENSOR(2,1,3,2) = DMATRIX(4,5)/COE2
      TENSOR(2,1,1,3) = DMATRIX(4,6)/COE2
      TENSOR(2,1,3,1) = DMATRIX(4,6)/COE2
C
      TENSOR(3,2,1,1) = DMATRIX(5,1)/COE1
      TENSOR(3,2,2,2) = DMATRIX(5,2)/COE1
      TENSOR(3,2,3,3) = DMATRIX(5,3)/COE1
      TENSOR(3,2,1,2) = DMATRIX(5,4)/COE2
      TENSOR(3,2,2,1) = DMATRIX(5,4)/COE2
      TENSOR(3,2,2,3) = DMATRIX(5,5)/COE2
      TENSOR(3,2,3,2) = DMATRIX(5,5)/COE2
      TENSOR(3,2,1,3) = DMATRIX(5,6)/COE2
      TENSOR(3,2,3,1) = DMATRIX(5,6)/COE2
C
      TENSOR(3,1,1,1) = DMATRIX(6,1)/COE1
      TENSOR(3,1,2,2) = DMATRIX(6,2)/COE1
      TENSOR(3,1,3,3) = DMATRIX(6,3)/COE1
      TENSOR(3,1,1,2) = DMATRIX(6,4)/COE2
      TENSOR(3,1,2,1) = DMATRIX(6,4)/COE2
      TENSOR(3,1,2,3) = DMATRIX(6,5)/COE2
      TENSOR(3,1,3,2) = DMATRIX(6,5)/COE2
      TENSOR(3,1,1,3) = DMATRIX(6,6)/COE2
      TENSOR(3,1,3,1) = DMATRIX(6,6)/COE2
C
C
      RETURN
      END
C
      SUBROUTINE IDENTITYMAT2(N,EYE2)
C========================================================================
C                                                                       =
C                     RANK-N ITENTITY MATRIX                            =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION EYE2(N,N)
      DATA ZERO,ONE /0.0D0,1.0D0/
C         
      DO I=1,N
          DO J=1,N
              EYE2(I,J)=ZERO
              IF (I.EQ.J) EYE2(I,J)=ONE
          END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE IDENTITYMAT4(N,EYE4)
C========================================================================
C                                                                       =
C                     RANK-N ITENTITY MATRIX                            =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION EYE4(N,N,N,N),EYE2(N,N)
      DATA ZERO,ONE /0.0D0,1.0D0/
C         
      CALL IDENTITYMAT2(3,EYE2)
      DO I=1,N
          DO J=1,N
              DO K=1,N
                  DO L=1,N
                      EYE4(I,J,K,L)=0.5*(EYE2(I,K)*EYE2(J,L)
     2                +EYE2(I,L)*EYE2(J,K))
                  END DO
              END DO
          END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aij_Bj(N,A,B,C)
C========================================================================
C                                                                       =
C                              Aij_Bj                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(N,N),B(N),C(N)
      DATA ZERO /0.0D0/
      DO I=1,N
        C(I)=ZERO
        DO J=1,N
          C(I)=C(I)+A(I,J)*B(J)
        END DO
      END DO
C
      RETURN
      END
      SUBROUTINE Aijkl_Bklmn(NN,A,B,C)
C========================================================================
C                                                                       =
C                              Aijkl_Bkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(NN,NN,NN,NN),B(NN,NN,NN,NN),C(NN,NN,NN,NN)
      DATA ZERO /0.0D0/
      DO I=1,NN
        DO J=1,NN
          DO K=1,NN
            DO L=1,NN
              C(I,J,K,L)=ZERO
              DO M=1,NN
                DO N=1,NN
                  C(I,J,K,L)=C(I,J,K,L)+A(I,J,M,N)*B(M,N,K,L)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
!C
      SUBROUTINE Aijkl_Bkl(N,A,B,C)
C========================================================================
C                                                                       =
C                              Aijkl_Bkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(N,N,N,N),B(N,N),C(N,N)
      DATA ZERO /0.0D0/
C
      DO I=1,N
        DO J=1,N
          C(I,J)=ZERO
          DO K=1,N
            DO L=1,N
              C(I,J)=C(I,J)+A(I,J,K,L)*B(K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
!C
      SUBROUTINE Aij_Bkl(N,A,B,C)
C====================================================================================
C                                                                                   =
C                                   Aij_Bkl                                         =
C                                                                                   =
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(N,N),B(N,N),C(N,N,N,N)
C
      DO I=1,N
          DO J=1,N
              DO K=1,N
                  DO L=1,N
                      C(I,J,K,L)=A(I,J)*B(K,L)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE Aik_Bjl(N,A,B,C)
C====================================================================================
C                                                                                   =
C                                   Aik_Bjl                                         =
C                                                                                   =
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(N,N),B(N,N),C(N,N,N,N)
C
      DO I=1,N
          DO J=1,N
              DO K=1,N
                  DO L=1,N
                      C(I,J,K,L)=A(I,K)*B(J,L)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE Ai_Bj(N,A,B,C)
C====================================================================================
C                                                                                   =
C                                   Ai_Bj                                         =
C                                                                                   =
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(N),B(N),C(N,N)
C
      DO I=1,N
        DO J=1,N
          C(I,J)=A(I)*B(J)
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Ai_Bj_Ck_Dl(N,A,B,C,D,E)
C====================================================================================
C                                                                                   =
C                                   Ai_Bj_Ck_Dl                                     =
C                                                                                   =
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(N),B(N),C(N),D(N),E(N,N,N,N)
C
      DO I=1,N
        DO J=1,N
          DO K=1,N
              DO L=1,N
                  E(I,J,K,L)=A(I)*B(J)*C(K)*D(L)
              END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C
C
      SUBROUTINE Aik_Bkj(N,A,B,C)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(N,N),B(N,N),C(N,N)
C
C
      DO I=1,N
        DO J=1,N
          C(I,J)=0.0D0
          DO K=1,N
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aij_Bij(N,A,B,C)
CC====================================================================================================
C
C                          Aij_Bij=A:B=C  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(N,N),B(N,N)
C
C     
      C=0.D0
      DO I=1,N
        DO J=1,N
          C=C+A(I,J)*B(I,J)
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aijmn_Bmnkl(A,B,C)
CC====================================================================================================
C
C                          Aijmnk_Bmnkl=Cijkl	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C(I,J,K,L)=0.0D0
              DO M=1,3
                DO N=1,3
            C(I,J,K,L)=C(I,J,K,L)+A(I,J,M,N)*B(M,N,K,L)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE MAT2_MAT1(NTENS,TENSOR,VECTOR,FACT)
C
C ====================================================================
C
C =================== MAT2_MAT1: TENSOR TO VECTOR=====================
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      VECTOR=0.0D0
C
      VECTOR( 1 ) = TENSOR(1 , 1)
      VECTOR( 2 ) = TENSOR(2 , 2)
      VECTOR( 3 ) = TENSOR(3 , 3)
      VECTOR( 4 ) = TENSOR(1 , 2)*FACT
      IF (NTENS>4) THEN
          VECTOR( 5 ) = TENSOR(1 , 3)*FACT
          VECTOR( 6 ) = TENSOR(2 , 3)*FACT
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE INVERSE(A,N,NP,AINV)
C========================================================================
C
C    CALCULATE THE SECOND ORDER TENSOR A'S INVERSE, AINV
C    A^{-1} = AINV    
C    this subroutine inverses a (n x n) A matrix
C	 following a Gauss-Jordan elimination process
C
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C  
      DIMENSION A(NP,NP),IPIV(NP),INDXR(NP),INDXC(NP),
     1 A0(NP,NP),AINV(NP,NP)
C
      DO J=1,N
        IPIV(J)=0
      END DO
C
C
C     storage of the original A matrix
C
      DO I=1,N
        DO J=1,N
          A0(I,J)=A(I,J)
        END DO
      END DO
C
C	find a pivot among the rows of A that have not already been reduced
C
      DO I=1,N
        BIG=0.0D0
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
                IF(IPIV(K).EQ.0)THEN
                  IF(ABS(A(J,K)).GE.BIG)THEN
                   BIG=ABS(A(J,K))
                   IROW=J
                   ICOL=K
                   PIV=A(J,K)
                  END IF
                ELSEIF(IPIV(K).GT.1)THEN
                  write (7,*) 'Singular Matrix'
              END IF
            END DO
          END IF
        END DO
C
      IPIV(ICOL)=IPIV(ICOL)+1
      INDXR(I)=IROW
      INDXC(I)=ICOL
C	  
C     interchange the rows to put the pivot on the diagonal
C
      IF(irow.ne.icol)THEN
        DO L=1,N
          DUM=A(IROW,L)
          A(IROW,L)=A(ICOL,L)
          A(ICOL,L)=DUM
        END DO
      END IF
C     reduction of the row of the pivot
C       
      IF(PIV.EQ.0) write (7,*) 'Singular Matrix2'
C       
      PIVIN=1./PIV          ! numerical stabilization
C
      A(ICOL,ICOL)=1.       ! numerical stabilization
        DO L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVIN
        END DO
C
C     reduction of the column of the pivot
C
        DO LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.    ! numerical stabilization
            DO L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
            END DO
          END IF
        END DO
      END DO
C
C
C     unscramble the columns to get A-1
C		
      DO J=N,1,-1   ! reverse DO loop
        DO K=1,N
          DUM=A(K,INDXR(J))
          A(K,INDXR(J))=A(K,INDXC(J))
          A(K,INDXC(J))=DUM
        END DO
      END DO
C
C	restitution process of A and Ainv
C
      DO I=1,N
        DO J=1,N
          AINV(I,J)=A(I,J)
        END DO
      END DO
C
      DO I=1,N
        DO J=1,N
          A(I,J)=A0(I,J)
        END DO
      END DO
C     
      RETURN
      END
C
      SUBROUTINE GEIGEN( N,A0,WR,WI,Z,FV1,IV1)
C=======================================================================
C
C      CALCULATE EIGENVALUES AND EIGENVECTORS OF A MATRIX
C    After https://ww2.odu.edu/~agodunov/computing/programs/others/Geigen.f      
C   Input
C	N - dimension
C	A0 - input array
C   Output
C     WR  - array of eigenvalues (real)
C	WI  - array of eigenvalues (imagin)
C	Z   - array of eigenvectors
C	FV1 - working array
C	IV1 - working array
C======================================================================
CCCCCCCOMMON/ KONTR/KTRW(80)
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A0(N,N)
      INTEGER*4 IV1(N)
      REAL*8 A(N,N),WR(N),WI(N),Z(N,N),FV1(N)
CCC
CCC   SUBSTITUTE A0 INTO ARRAY A TO KEEP INPUT INTACT
      DO I=1,N
          DO J=1,N
              A(I,J)=A0(I,J)
          END DO
      END DO
CCC   SUBROUTINES FROM EISPECK-PACKAGE
      CALL BALANC(N,N,A, IS1,IS2,FV1)
      CALL ELMHES(N,N,IS1,IS2,A,IV1)
      CALL GEIELTRAN(N,N,IS1,IS2,A,IV1,Z)
      CALL HQR2(N,N,IS1,IS2,A,WR,WI,Z,IERR)
CG    IF(IERR.NE.0)  CALL WRITEI('IERR',1,IERR,1)
      CALL BALBAK(N,N,IS1,IS2,FV1,N,Z)
C
CCC---
CG    IF(KTRW(4).EQ.1)   CALL WRITED('WR  ',1,WR,N)
CG    IF(KTRW(4).EQ.1)   CALL WRITED('WI  ',1,WI,N)
      RETURN
      END
C
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
C===============================================================================
C===============================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL*8 A(NM,N),SCALE(N)
      REAL*8 C,F,G,R,S,B2,RADIX
      REAL*8 DABS
      LOGICAL NOCONV
C
C     ***** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C          THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION
C     ******
      RADIX=16
C
      B2=RADIX*RADIX
      K=1
      L=N
      GO TO 100
C     ***** IN-LINE PROCEDURE FOR ROW AND
C          COLUMN EXCHANGE*****
20    SCALE(M)=J
      IF(J.EQ.M) GO TO 50
C
       DO 30 I=1,L
      F=A(I,J)
      A(I,J)=A(I,M)
      A(I,M)=F
30    CONTINUE
C
       DO  40 I=K,N
      F=A(J,I)
      A(J,I)=A(M,I)
      A(M,I)=F
40    CONTINUE
C
50    GO TO (80,130),IEXC
C     *****SEARCH FOR ROWS IS_EFF0LATING AN EIGENVALUE
C          AND PUSH THEM DOWN*****
80    IF(L.EQ.1) GO TO 280
      L=L-1
C     *****FOR J=L STEP -1 UNTIL 1 DO -- *****
100   DO 120 JJ=1,L
      J=L+1-JJ
C
       DO 110 I=1,L
      IF(I.EQ.J) GO TO 110
      IF(A(J,I).NE.0.0) GO TO 120
110   CONTINUE
C
      M=L
      IEXC=1
      GO TO 20
120   CONTINUE
C
      GO TO 140
C     *****SEARCH FOR COLUMS IS_EFF0LATING AN EIGENVALUE
C          AND PUSH THEM LEFT*****
130   K=K+1
C
140   DO 170 J=K,L
C
       DO 150 I=K,L
      IF(I.EQ.J) GO TO 150
      IF(A(I,J).NE.0.0) GO TO 170
150   CONTINUE
C
      M=K
      IEXC=2
      GO TO 20
170   CONTINUE
C     *****NOW BALANCE THE SUBMATRIX IN ROWS K TO L*****
       DO 180 I=K,L
180   SCALE(I)=1.0
C     *****ITERATIVE LOOP FOR NORM REDUCTION*****
190   NOCONV=.FALSE.
C
       DO 270 I=K,L
      C=0.0
      R=0.0
C
       DO 200 J=K,L
      IF(J.EQ.I) GO TO 200
      C=C+DABS(A(J,I))
      R=R+DABS(A(I,J))
200   CONTINUE
C
      G=R/RADIX
      F=1.0
      S=C+R
210   IF(C.GE.G) GO TO 220
      F=F*RADIX
      C=C*B2
      GO TO 210
220   G=R*RADIX
230   IF(C.LT.G) GO TO 240
      F=F/RADIX
      C=C/B2
      GO TO 230
C     *****NOW BALANCE*****
240   IF((C+R)/F.GE.0.95*S) GO TO 270
      G=1.0/F
      SCALE(I)=SCALE(I)*F
      NOCONV=.TRUE.
C
       DO 250 J=K,N
250   A(I,J)=A(I,J)*G
C
       DO 260 J=1,L
260   A(J,I)=A(J,I)*F
C
270   CONTINUE
C
      IF(NOCONV) GO TO 190
C
280   LOW=K
      IGH=L
      RETURN
      END
C
C   02/03/76            MEMBER NAME  ELMHES   (QUELLE)      FORTRAN
      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
C=================================================================================
C=================================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,NM1,MP1
      REAL*8 A(NM,N)
      REAL*8 X,Y
      REAL*8 DABS
      INTEGER INT(IGH)
C
      LA=IGH-1
      KP1=LOW+1
      IF(LA.LT.KP1) GO TO 200
C
       DO 180 M=KP1,LA
      MM1=M-1
      X=0.0
      I=M
C
       DO 100 J=M,IGH
      IF(DABS(A(J,MM1)).LE.DABS(X)) GO TO 100
      X=A(J,MM1)
      I=J
100   CONTINUE
C
      INT(M)=I
      IF(I.EQ.M) GO TO 130
C     *****INTERCHANGE ROWS AND COLUMNS OF A*****
       DO 110 J=MM1,N
      Y=A(I,J)
      A(I,J)=A(M,J)
      A(M,J)=Y
110   CONTINUE
C
       DO 120 J=1,IGH
      Y=A(J,I)
      A(J,I)=A(J,M)
      A(J,M)=Y
120   CONTINUE
C     *****END INTERCHANGE*****
130   IF(X.EQ.0.0) GO TO 180
      MP1=M+1
C
       DO 160 I=MP1,IGH
      Y=A(I,MM1)
      IF(Y.EQ.0.0) GO TO 160
      Y=Y/X
      A(I,MM1)=Y
C
       DO 140 J=M,N
140   A(I,J)=A(I,J)-Y*A(M,J)
C
       DO 150 J=1,IGH
150   A(J,M)=A(J,M)+Y*A(J,I)
C
160   CONTINUE
C
180   CONTINUE
C
200   RETURN
      END
C   02/03/76            MEMBER NAME  ELTRAN   (QUELLE)      FORTRAN
      SUBROUTINE GEIELTRAN(NM,N,LOW,IGH,A,INT,Z)
C======================================================================================
C-======================================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,N,KL,NM,MP,NM1,IGH,LOW,MP1
      REAL*8 A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
C
C     *****INITIALIZE Z TO IDENTITY MATRIX*****
       DO 80 I=1,N
C
       DO 60 J=1,N
60    Z(I,J)=0.0
C
      Z(I,I)=1.0
80    CONTINUE
C
      KL=IGH-LOW-1
      IF(KL.LT.1) GO TO 200
C     *****FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- *****
       DO 140 MM=1,KL
      MP=IGH-MM
      MP1=MP+1
C
       DO 100 I=MP1,IGH
100   Z(I,MP)=A(I,MP-1)
C
      I=INT(MP)
      IF(I.EQ.MP) GO TO 140
C
       DO 130 J=MP,IGH
      Z(MP,J)=Z(I,J)
      Z(I,J)=0.0
130   CONTINUE
C
      Z(I,MP)=1.0
140   CONTINUE
C
200   RETURN
      END
C   02/03/76            MEMBER NAME  HQR2     (QUELLE)      FORTRAN
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C=====================================================================================
C=====================================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITS,LOW,MP2,ENM2,IERR
      REAL*8 H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL*8 P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
      REAL*8 DSQRT,DABS,DSIGN
      INTEGER MIN0
      LOGICAL NOTLAS
      COMPLEX*16 Z3
      COMPLEX*16 DCMPLX
      REAL*8 T3(2)
      EQUIVALENCE(Z3,T3(1))
C
C     *****MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C          THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C     *****
      MACHEP= 1.0D-16
C
      IERR=0
C     *****STORE ROOTS IS_EFF0LATED BY BALANC*****
       DO 50 I=1,N
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 50
      WR(I)=H(I,I)
      WI(I)=0.0
50    CONTINUE
C
      EN=IGH
      T=0.0
C     *****SEARCH FOR NEXT EIGENVALUES*****
60    IF(EN.LT.LOW) GO TO 340
      ITS=0
      NA=EN-1
      ENM2=NA-1
C     *****LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C          FOR L=EN STEP -1 UNTIL LOW DO -- *****
70    DO 80 LL=LOW,EN
      L=EN+LOW-LL
      IF(L.EQ.LOW) GO TO 100
      IF(DABS(H(L,L-1)).LE.MACHEP*(DABS(H(L-1,L-1))
     X    +DABS(H(L,L)))) GO TO 100
80    CONTINUE
C     *****FORM SHIFT*****
100   X=H(EN,EN)
      IF(L.EQ.EN) GO TO 270
      Y=H(NA,NA)
      W=H(EN,NA)*H(NA,EN)
      IF(L.EQ.NA) GO TO 280
      IF(ITS.EQ.30) GO TO 1000
      IF(ITS.NE.10.AND.ITS.NE.20) GO TO 130
C     *****FORM EXCEPTIONAL SHIFT*****
      T=T+X
C
       DO 120 I=LOW,EN
120   H(I,I)=H(I,I)-X
C
      S=DABS(H(EN,NA))+DABS(H(NA,ENM2))
      X=0.75*S
      Y=X
      W=-0.4375*S*S
130   ITS=ITS+1
C     LOOK FOR TWO CONSECUTIVE SMALL
C          SUB-DIAGONAL ELEMENTS.
C          FOR M=EN-2 STEP -1 UNTIL L DO -- *****
       DO 140 MM=L,ENM2
      M=ENM2+L-MM
      ZZ=H(M,M)
      R=X-ZZ
      S=Y-ZZ
      P=(R*S-W)/H(M+1,M)+H(M,M+1)
      Q=H(M+1,M+1)-ZZ-R-S
      R=H(M+2,M+1)
      S=DABS(P)+DABS(Q)+DABS(R)
      P=P/S
      Q=Q/S
      R=R/S
      IF(M.EQ.L) GO TO 150
      IF(DABS(H(M,M-1))*(DABS(Q)+DABS(R)).LE.MACHEP*DABS(P)
     X  *(DABS(H(M-1,M-1))+DABS(ZZ)+DABS(H(M+1,M+1)))) GO TO 150
140   CONTINUE
C
150   MP2=M+2
C
       DO 160 I=MP2,EN
      H(I,I-2)=0.0
      IF(I.EQ.MP2) GO TO 160
      H(I,I-3)=0.0
160   CONTINUE
C     *****DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C          COLUMNS M TO EN*****
       DO 260 K=M,NA
      NOTLAS=K.NE.NA
      IF(K.EQ.M) GO TO 170
      P=H(K,K-1)
      Q=H(K+1,K-1)
      R=0.0
      IF(NOTLAS) R=H(K+2,K-1)
      X=DABS(P)+DABS(Q)+DABS(R)
      IF(X.EQ.0.0) GO TO 260
      P=P/X
      Q=Q/X
      R=R/X
170   S=DSIGN(DSQRT(P*P+Q*Q+R*R),P)
      IF(K.EQ.M) GO TO 180
      H(K,K-1)=-S*X
      GO TO 190
180   IF(L.NE.M) H(K,K-1)=-H(K,K-1)
190   P=P+S
      X=P/S
      Y=Q/S
      ZZ=R/S
      Q=Q/P
      R=R/P
C     *****ROW MODIFICATION*****
       DO 210 J=K,N
      P=H(K,J)+Q*H(K+1,J)
      IF(.NOT.NOTLAS) GO TO 200
      P=P+R*H(K+2,J)
      H(K+2,J)=H(K+2,J)-P*ZZ
200   H(K+1,J)=H(K+1,J)-P*Y
      H(K,J)=H(K,J)-P*X
210   CONTINUE
C
      J=MIN0(EN,K+3)
C     *****COLUMN MODIFICATION*****
       DO 230 I=1,J
      P=X*H(I,K)+Y*H(I,K+1)
      IF(.NOT.NOTLAS) GO TO 220
      P=P+ZZ*H(I,K+2)
      H(I,K+2)=H(I,K+2)-P*R
220   H(I,K+1)=H(I,K+1)-P*Q
      H(I,K)=H(I,K)-P
230   CONTINUE
C     *****ACCUMULATE TRANSFORMATIONS*****
       DO 250 I=LOW,IGH
      P=X*Z(I,K)+Y*Z(I,K+1)
      IF(.NOT.NOTLAS) GO TO 240
      P=P+ZZ*Z(I,K+2)
      Z(I,K+2)=Z(I,K+2)-P*R
240   Z(I,K+1)=Z(I,K+1)-P*Q
      Z(I,K)=Z(I,K)-P
250   CONTINUE
C
260   CONTINUE
C
      GO TO 70
C     *****ONE ROOT FOUND*****
270   H(EN,EN)=X+T
      WR(EN)=H(EN,EN)
      WI(EN)=0.0
      EN=NA
      GO TO 60
C     *****TWO ROOTS FOUND*****
280   P=(Y-X)/2.0
      Q=P*P+W
      ZZ=DSQRT(DABS(Q))
      H(EN,EN)=X+T
      X=H(EN,EN)
      H(NA,NA)=Y+T
      IF(Q.LT.0.0) GO TO 320
C     *****REAL PAIR*****
      ZZ=P+DSIGN(ZZ,P)
      WR(NA)=X+ZZ
      WR(EN)=WR(NA)
      IF(ZZ.NE.0.0) WR(EN)=X-W/ZZ
      WI(NA)=0.0
      WI(EN)=0.0
      X=H(EN,NA)
      R=DSQRT(X*X+ZZ*ZZ)
      P=X/R
      Q=ZZ/R
C     *****ROW MODIFICATION*****
       DO 290 J=NA,N
      ZZ=H(NA,J)
      H(NA,J)=Q*ZZ+P*H(EN,J)
      H(EN,J)=Q*H(EN,J)-P*ZZ
290   CONTINUE
C     *****COLUMN MODIFICATION*****
       DO 300 I=1,EN
      ZZ=H(I,NA)
      H(I,NA)=Q*ZZ+P*H(I,EN)
      H(I,EN)=Q*H(I,EN)-P*ZZ
300   CONTINUE
C     *****ACCUMULATE TRANSFORMATIONS*****
       DO 310 I=LOW,IGH
      ZZ=Z(I,NA)
      Z(I,NA)=Q*ZZ+P*Z(I,EN)
      Z(I,EN)=Q*Z(I,EN)-P*ZZ
310   CONTINUE
C
      GO TO 330
C     *****COMPLEX PAIR*****
320   WR(NA)=X+P
      WR(EN)=X+P
      WI(NA)=ZZ
      WI(EN)=-ZZ
330   EN=ENM2
      GO TO 60
C     *****ALL ROOTS FOUND.BACKSUBSTITUTE TO FIND
C          VECTORS OF UPPER TRIANGULAR FORM*****
340   NORM=0.0
      K=1
C
       DO 360 I=1,N
C
       DO 350 J=K,N
350   NORM=NORM+DABS(H(I,J))
C
      K=I
360   CONTINUE
C
      IF(NORM.EQ.0.0) GO TO 1001
C     *****FOR EN=N STEP -1 UNTIL 1 DO -- *****
       DO 800 NN=1,N
      EN=N+1-NN
      P=WR(EN)
      Q=WI(EN)
      NA=EN-1
      IF(Q) 710,600,800
C     *****REAL VECTOR*****
600   M=EN
      H(EN,EN)=1.0
      IF(NA.EQ.0) GO TO 800
C     *****FOR I=EN-1 STEP -1 UNTIL 1 DO --*****
       DO 700 II=1,NA
      I=EN-II
      W=H(I,I)-P
      R=H(I,EN)
      IF(M.GT.NA) GO TO 620
C
       DO 610 J=M,NA
610   R=R+H(I,J)*H(J,EN)
C
620   IF(WI(I).GE.0.0) GO TO 630
      ZZ=W
      S=R
      GO TO 700
630   M=I
      IF(WI(I).NE.0.0) GO TO 640
      T=W
      IF(W.EQ.0.0) T=MACHEP*NORM
      H(I,EN)=-R/T
      GO TO 700
C************* SOLVE REAL EQUATIONS *******************
640   X=H(I,I+1)
      Y=H(I+1,I)
      Q=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
      T=(X*S-ZZ*R)/Q
      H(I,EN)=T
      IF(DABS(X).LE.DABS(ZZ)) GOTO 650
      H(I+1,EN)=(-R-W*T)/X
      GOTO 700
650   H(I+1,EN)=(-S-Y*T)/ZZ
700   CONTINUE
C     *****END REAL VECTOR*****
      GO TO 800
C     *****COMPLEX VECTOR*****
710   M=NA
C     *****LAST VECTOR COMPONENT CHOSEN IMAGINARY S_EFF0 THAT
C          EIGENVECTOR MATRIX IS TRIANGULAR*****
      IF(DABS(H(EN,NA)).LE.DABS(H(NA,EN))) GO TO 720
      H(NA,NA)=Q/H(EN,NA)
      H(NA,EN)=-(H(EN,EN)-P)/H(EN,NA)
      GO TO 730
720   Z3=DCMPLX(0.D+0,-H(NA,EN))/DCMPLX(H(NA,NA)-P,Q)
      H(NA,NA)=T3(1)
      H(NA,EN)=T3(2)
730   H(EN,NA)=0.0
      H(EN,EN)=1.0
      ENM2=NA-1
      IF(ENM2.EQ.0) GO TO 800
C
       DO 790 II=1,ENM2
      I=NA-II
      W=H(I,I)-P
      RA=0.0
      SA=H(I,EN)
C
       DO 760 J=M,NA
      RA=RA+H(I,J)*H(J,NA)
      SA=SA+H(I,J)*H(J,EN)
760   CONTINUE
C
      IF(WI(I).GE.0.0) GO TO 770
      ZZ=W
      R=RA
      S=SA
      GO TO 790
770   M=I
      IF(WI(I).NE.0.0) GO TO 780
      Z3=DCMPLX(-RA,-SA)/DCMPLX(W,Q)
      H(I,NA)=T3(1)
      H(I,EN)=T3(2)
      GO TO 790
C     *****S_EFF0LVE COMPLEX EQUATIONS*****
780   X=H(I,I+1)
      Y=H(I+1,I)
      VR=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
      VI=(WR(I)-P)*2.0*Q
      IF(VR.EQ.0.0.AND.VI.EQ.0.0) VR=MACHEP*NORM
     X  *(DABS(W)+DABS(Q)+DABS(X)+DABS(Y)+DABS(ZZ))
      Z3=DCMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/DCMPLX(VR,VI)
      H(I,NA)=T3(1)
      H(I,EN)=T3(2)
      IF(DABS(X).LE.DABS(ZZ)+DABS(Q)) GO TO 785
      H(I+1,NA)=(-RA-W*H(I,NA)+Q*H(I,EN))/X
      H(I+1,EN)=(-SA-W*H(I,EN)-Q*H(I,NA))/X
      GO TO 790
785   Z3=DCMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN))/DCMPLX(ZZ,Q)
      H(I+1,NA)=T3(1)
      H(I+1,EN)=T3(2)
790   CONTINUE
C     *****END COMPLEX VECTOR*****
800   CONTINUE
C     *****END BACK SUBSTITUTION.
C          VECTORS OF IS_EFF0LATED ROOTS*****
       DO 840 I=1,N
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 840
C
       DO 820 J=1,N
820   Z(I,J)=H(I,J)
C
840   CONTINUE
C     *****MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C          VECTORS OF ORIGINAL FULL MATRIX.
C          FOR J=N STEP -1 UNTIL LOW DO --*****
       DO 880 JJ=LOW,N
      J=N+LOW-JJ
      M=MIN0(J,IGH)
C
       DO 880 I=LOW,IGH
      ZZ=0.0
C
       DO 860 K=LOW,M
860   ZZ=ZZ+Z(I,K)*H(K,J)
C
      Z(I,J)=ZZ
880   CONTINUE
C
      GO TO 1001
C     *****SET ERROR -- NO CONVERGENCE TO AN
C          EIGENVALUE AFTER 30 ITERATIONS*****
1000  IERR=EN
1001  RETURN
      END
C   02/03/76            MEMBER NAME  BALBAK   (QUELLE)      FORTRAN
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
C========================================================================
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL*8 SCALE(N),Z(NM,M)
      REAL*8 S,SUM
C
      IF(IGH.EQ.LOW) GO TO 120
C
       DO 110 I=LOW,IGH
      S=SCALE(I)
C     *****LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C          IF THE FOREGOING STATEMENT IS REPLACED BY
C          S=1.0/SCALE(I).*****
       DO 100 J=1,M
100   Z(I,J)=Z(I,J)*S
C
110   CONTINUE
C     *****FOR I=LOW-1 STEP -1 UNTIL 1,
C          IGH+1 STEP 1 UNTIL N DO -- *****
120   DO 140 II=1,N
      I=II
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 140
      IF(I.LT.LOW) I=LOW-II
      K=SCALE(I)
      IF(K.EQ.I) GO TO 140
C
       DO 130 J=1,M
      S=Z(I,J)
      Z(I,J)=Z(K,J)
      Z(K,J)=S
130   CONTINUE
C
140   CONTINUE
C
       DO  170 J=1,M
      SUM=0.0
       DO  150 I=1,N
      SUM=SUM+Z(I,J)*Z(I,J)
150   CONTINUE
      SUM=DSQRT(SUM)
       DO  160 I=1,N
      Z(I,J)=Z(I,J)/SUM
160   CONTINUE
170   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE TRACE(N,Aij,Aii)
C==================================================================
C
C             FIND TRACE OF A RANK-N MATRIX
C
C===================================================================
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION Aij(N,N)
      Aii=0.0D0
      DO i=1,N
          Aii=Aii+Aij(i,i)
      END DO
      RETURN
      END
C
C
      SUBROUTINE POS_PROJ(EPS,PPLUS,EPSPLUS)
C==================================================================
C
C             FIND POSITIVE PROJECTION OF STRAIN TENSOR
C
C===================================================================
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION EPS(3,3),IV1(3),EPSPLUS(3,3),EIGVEC_i(3),EIGVEC_j(3),
     2 EIGVEC_ij(3,3),EPS_PLUS(3,3),QPLUS(3,3),PPLUS(3,3,3,3)
      REAL*8 FV1(3),EIGVALR(3),EIGVALI(3),EIGVEC(3,3)
      DO I=1,3
          FV1(I)=0.
          IV1(I)=0.
          EIGVALR(I)=0.
          EIGVALI(I)=0.
          DO J=1,3
              EPSPLUS(I,J)=0.
              QPLUS(I,J)=0.
              EIGVEC(I,J)=0.
          END DO
      END DO
      CALL GEIGEN(3,EPS,EIGVALR,EIGVALI,EIGVEC,FV1,IV1)
      DO K=1,3
          HEAVISIDE=0
          IF (EIGVALR(K).GE.1.E-7) THEN
              HEAVISIDE=1
          ENDIF
          EIGVEC_i=EIGVEC(:,K)
          CALL Ai_Bj(3,EIGVEC_i,EIGVEC_i,EIGVEC_ij)
          EPSPLUS=EPSPLUS+HEAVISIDE*EIGVALR(K)*EIGVEC_ij
          QPLUS=QPLUS+HEAVISIDE*EIGVEC_ij
          CALL Aik_Bjl(3,QPLUS,QPLUS,PPLUS)
      END DO
      RETURN
      END 
C     
C
C
C
      SUBROUTINE D_FCN(EPS,D,FD)
C==================================================================
C
C             DAMAGE CRITERIA
C
C===================================================================
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION EPS(3,3),EPSPLUS(3,3),PPLUS(3,3,3,3)
C
      COMMON/PARAMETERS/ELAMBDA0,EMU0,CKAPPA,C0D,C1D,C2D
C
      CALL POS_PROJ(EPS,PPLUS,EPSPLUS)
      CALL Aij_Bij(3,EPSPLUS,EPSPLUS,EPSPLUS2)   
      FD=C1D-C1D*EXP((C0D-SQRT(EPSPLUS2))/C2D)-D
      IF (FD.GT.0) THEN
          D=C1D-C1D*EXP((C0D-SQRT(EPSPLUS2))/C2D)
      ENDIF
      RETURN
      END
C      
C
      SUBROUTINE CAL_C(NTENS,D,R,C,CNN)
C==================================================================
C
C             UPDATE STRESS & STIFFNESS TENSOR
C
C===================================================================
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION EYE2(3,3),
     3 C0(3,3,3,3),CD(3,3,3,3),CD99(9,9),CNN(NTENS,NTENS),
     4 C(3,3,3,3)
C
      COMMON/PARAMETERS/ELAMBDA0,EMU0,CKAPPA,C0D,C1D,C2D
C
      CALL IDENTITYMAT2(3,EYE2)
C     STIFFNESS TENSOR 
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
                      C0(I,J,K,L)=ELAMBDA0*EYE2(K,L)*EYE2(I,J)+
     2                EMU0*(EYE2(I,K)*EYE2(J,L)+EYE2(I,L)*EYE2(J,K))
                  END DO
              END DO
          END DO
      END DO
      C=(1-(1-R)*D)*C0   
C
      CALL MAT4_MAT2(NTENS,C,CNN,1)
      RETURN
      END
C
C