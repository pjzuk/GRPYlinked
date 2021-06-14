C February 10, 2015

C LAPACK necessary for SUBROUTINE EIGEN(A,D,NN)
C SUBROUTINE ZINIT returns Z ( and NOT Z^{-1} ).

C EVALUATES THE 11NNx11NN ROTNE-PRAGER MOBILITY MATRIX.
C DISPLAYS BLOKS:
C PP  6NNx6NN
C PQ  6NNx5NN
C QQ  5NNx5NN

C WORKS FOR ARBITRARY LL BECAUSE
C      H(B,L,M)=3*(L-1)*(L+1) + B*(2*L+1) + L+1+M
C IS USED AND NOT
C      H(B,L,M)=B*LLF + L*(L+1)+M .

*****************************************

      SUBROUTINE HYDRO(A,C,CONF,RADII)
      USE SIZE       ! NN,NN6,LLG
      USE MULT       ! LL
      USE LATTICE_SKEW
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A(3*NN,3*NN),C(3,NN)
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN),ARR(6*NN,6*NN)
      REAL*8 CONF(3,NN),RADII(NN),SQRT2
      INTEGER I,J

      SQRT2=SQRT(2.D0)

      A=0.D0
      C=0.D0

      CALL HYDRO_RP(APP,APQ,AQQ,CONF,RADII)

      ARR = APP
      CALL INVFRI_TO_FRI(ARR,APQ,AQQ,NN)
      CALL FRI_TO_MOB_RED(APP,APQ,AQQ,NN)

      DO I=1,NN
       DO J=1,NN
        A(3*(I-1)+1:3*I,3*(J-1)+1:3*J) = 
     * APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+1:6*(J-1)+3)
       ENDDO
      ENDDO

      DO I=1,NN
       DO J=1,NN
        C(1:3,I) = C(1:3,I) -
     * APQ(6*(I-1)+1:6*(I-1)+3,5*(J-1)+4)/SQRT2 
       ENDDO
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE HYDRO_RP(APP,APQ,AQQ,CONF,RADII)
      USE SIZE
      IMPLICIT NONE
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 CONF(3,NN),RADII(NN)

      CALL ROTNE_PRAGER_YAMAKAWA(APP,APQ,AQQ,CONF,RADII,NN)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_YAMAKAWA(APP,APQ,AQQ,CONF,RADII,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 CONF(3,NN),RADII(NN)
      REAL*8 A1(3,3),B1(5,5),C1(3,5),R(3),DIST
      REAL*8 aaI,aaJ            ! radius
      INTEGER I,J
      LOGICAL :: INIT=.TRUE.
      SAVE INIT
      
      IF(INIT) THEN
       INIT=.FALSE.
       CALL INIT_U
       CALL INIT_EPS
       CALL INIT_Y2
      ENDIF

      APP=0.D0
      APQ=0.D0
      AQQ=0.D0

C PARTICLES LOOP ********************************************

      DO I=1,NN
        aaI=RADII(I)
       DO J=1,3
        APP(6*(I-1)+J,6*(I-1)+J) = 1.D0/(6.D0*aaI)
        APP(6*(I-1)+3+J,6*(I-1)+3+J) = 1.D0/(8.D0*aaI**3)
       ENDDO
       DO J=1,5
        AQQ(5*(I-1)+J,5*(I-1)+J) = 3.D0/(20.D0*aaI**3)
       ENDDO
      ENDDO

      DO I=1,NN-1
      DO J=I+1,NN       
       
       R=CONF(1:3,I)-CONF(1:3,J)
       DIST = SQRT(SUM(R**2))       

       aaI=RADII(I)
       aaJ=RADII(J)
       
       IF(DIST>aaI+aaJ) THEN
      
       CALL ROTNE_PRAGER_TT_IJ(A1,R,aaI,aaJ)
       APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+1:6*(J-1)+3)=A1
       APP(6*(J-1)+1:6*(J-1)+3,6*(I-1)+1:6*(I-1)+3)=A1
      
       CALL ROTNE_PRAGER_RR_IJ(A1,R)
       APP(6*(I-1)+4:6*(I-1)+6,6*(J-1)+4:6*(J-1)+6)=A1
       APP(6*(J-1)+4:6*(J-1)+6,6*(I-1)+4:6*(I-1)+6)=A1
      
       CALL ROTNE_PRAGER_RT_IJ(A1,R)
       APP(6*(I-1)+4:6*(I-1)+6,6*(J-1)+1:6*(J-1)+3)=A1  ! RT IJ
       APP(6*(J-1)+4:6*(J-1)+6,6*(I-1)+1:6*(I-1)+3)=-A1 ! RT JI
       APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+4:6*(J-1)+6)=A1  ! TR IJ
       APP(6*(J-1)+1:6*(J-1)+3,6*(I-1)+4:6*(I-1)+6)=-A1 ! TR JI
       
       CALL ROTNE_PRAGER_TD_IJ_Y2(C1,R,aaI,aaJ)
       APQ(6*(I-1)+1:6*(I-1)+3,5*(J-1)+1:5*(J-1)+5)=-C1
       CALL ROTNE_PRAGER_TD_IJ_Y2(C1,-R,aaJ,aaI)
       APQ(6*(J-1)+1:6*(J-1)+3,5*(I-1)+1:5*(I-1)+5)=-C1
       
       CALL ROTNE_PRAGER_RD_IJ_Y2(C1,R,aaI,aaJ)
       APQ(6*(I-1)+4:6*(I-1)+6,5*(J-1)+1:5*(J-1)+5)=-C1
       CALL ROTNE_PRAGER_RD_IJ_Y2(C1,-R,aaJ,aaI)
       APQ(6*(J-1)+4:6*(J-1)+6,5*(I-1)+1:5*(I-1)+5)=-C1

       CALL ROTNE_PRAGER_DD_IJ_Y2(B1,R,aaI,aaJ)
       AQQ(5*(I-1)+1:5*(I-1)+5,5*(J-1)+1:5*(J-1)+5)=B1
       AQQ(5*(J-1)+1:5*(J-1)+5,5*(I-1)+1:5*(I-1)+5)=B1

C***************************************************************    
       ELSE
              
       CALL YAMAKAWA_TT_IJ(A1,R,aaI,aaJ)
       APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+1:6*(J-1)+3)=A1
       APP(6*(J-1)+1:6*(J-1)+3,6*(I-1)+1:6*(I-1)+3)=A1
       
       CALL YAMAKAWA_RR_IJ(A1,R,aaI,aaJ)
       APP(6*(I-1)+4:6*(I-1)+6,6*(J-1)+4:6*(J-1)+6)=A1
       APP(6*(J-1)+4:6*(J-1)+6,6*(I-1)+4:6*(I-1)+6)=A1
       
       CALL YAMAKAWA_RT_IJ(A1,R,aaI,aaJ)
       APP(6*(I-1)+4:6*(I-1)+6,6*(J-1)+1:6*(J-1)+3)=A1  ! RT IJ
       APP(6*(J-1)+1:6*(J-1)+3,6*(I-1)+4:6*(I-1)+6)=-A1 ! TR JI

       CALL YAMAKAWA_RT_IJ(A1,-R,aaJ,aaI)
       APP(6*(J-1)+4:6*(J-1)+6,6*(I-1)+1:6*(I-1)+3)=A1  ! RT JI
       APP(6*(I-1)+1:6*(I-1)+3,6*(J-1)+4:6*(J-1)+6)=-A1 ! TR IJ

       CALL YAMAKAWA_TD_IJ_Y2(C1,R,aaI,aaJ)
       APQ(6*(I-1)+1:6*(I-1)+3,5*(J-1)+1:5*(J-1)+5)=-C1
       CALL YAMAKAWA_TD_IJ_Y2(C1,-R,aaJ,aaI)
       APQ(6*(J-1)+1:6*(J-1)+3,5*(I-1)+1:5*(I-1)+5)=-C1
       
       CALL YAMAKAWA_RD_IJ_Y2(C1,R,aaI,aaJ)
       APQ(6*(I-1)+4:6*(I-1)+6,5*(J-1)+1:5*(J-1)+5)=-C1
       CALL YAMAKAWA_RD_IJ_Y2(C1,-R,aaJ,aaI)
       APQ(6*(J-1)+4:6*(J-1)+6,5*(I-1)+1:5*(I-1)+5)=-C1

       CALL YAMAKAWA_DD_IJ_Y2(B1,R,aaI,aaJ)
       AQQ(5*(I-1)+1:5*(I-1)+5,5*(J-1)+1:5*(J-1)+5)=B1
       AQQ(5*(J-1)+1:5*(J-1)+5,5*(I-1)+1:5*(I-1)+5)=B1
       
       ENDIF
       
      ENDDO
      ENDDO
      
      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE MAT_INV3(LI,LR,DET_LR)
      IMPLICIT NONE
      REAL*8 LI(3,3),LR(3,3),DET_LR
      
      det_LR=LR(1,1)*LR(2,2)*LR(3,3)+LR(1,2)*LR(2,3)*LR(3,1)+
     *        LR(1,3)*LR(2,1)*LR(3,2)-LR(1,3)*LR(2,2)*LR(3,1)-      
     *        LR(2,3)*LR(3,2)*LR(1,1)-LR(3,3)*LR(1,2)*LR(2,1)

      LI(1,1)=(1.D0/det_LR)*(LR(2,2)*LR(3,3)-LR(2,3)*LR(3,2))
      LI(1,2)=(1.D0/det_LR)*(LR(1,3)*LR(3,2)-LR(1,2)*LR(3,3))
      LI(1,3)=(1.D0/det_LR)*(LR(1,2)*LR(2,3)-LR(1,3)*LR(2,2))

      LI(2,1)=(1.D0/det_LR)*(LR(2,3)*LR(3,1)-LR(2,1)*LR(3,3))
      LI(2,2)=(1.D0/det_LR)*(LR(1,1)*LR(3,3)-LR(1,3)*LR(3,1))
      LI(2,3)=(1.D0/det_LR)*(LR(1,3)*LR(2,1)-LR(1,1)*LR(2,3))

      LI(3,1)=(1.D0/det_LR)*(LR(2,1)*LR(3,2)-LR(2,2)*LR(3,1))
      LI(3,2)=(1.D0/det_LR)*(LR(1,2)*LR(3,1)-LR(1,1)*LR(3,2))
      LI(3,3)=(1.D0/det_LR)*(LR(1,1)*LR(2,2)-LR(1,2)*LR(2,1))
      
      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

C NOT NORMALIZED
      SUBROUTINE ROTNE_PRAGER_TT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST3,RW(3),RR(3,3)
      REAL*8 aaI,aaJ           ! radius
      REAL*8 P18,P23a2
      PARAMETER(P18=1.D0/8.D0)
      P23a2=(aaI**2 + aaJ**2)/3.D0

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST3=DIST**3

      CALL CALC_RR(RR,RW)

      A1=P18*( (1.D0/DIST+P23a2/DIST3)*U + (1/DIST-3.D0*P23a2/DIST3)*RR)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

C NOT NORMALIZED
      SUBROUTINE YAMAKAWA_TT_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST2,DIST3,RW(3),RR(3,3)
      REAL*8 aaI,aaJ            ! radius
      REAL*8 M0TT          ! tt single sphere mobility 
      PARAMETER(M0TT=1.D0/6.D0)
      REAL*8 PRE,PU,PRR,aaIaaJ

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_RR(RR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=M0TT/aaJ*U
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        PRE=1.D0/(6*32*aaI*aaJ)
        PU=16*(aaI+aaJ)-(aaIaaJ**2+3*DIST2)**2/DIST3
        PRR=3*(aaIaaJ**2 - DIST2)**2/DIST3
        A1=PRE*(PU*U+PRR*RR)
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=M0TT/aaI*U
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        PRE=1.D0/(6*32*aaI*aaJ)
        PU=16*(aaI+aaJ)-(aaIaaJ**2+3*DIST2)**2/DIST3
        PRR=3*(aaIaaJ**2 - DIST2)**2/DIST3
        A1=PRE*(PU*U+PRR*RR)
       ENDIF       
      ENDIF

      RETURN
      END
      
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_RR_IJ(A1,R)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST3,RW(3),RR(3,3)
      REAL*8 P16
      PARAMETER(P16=1.D0/16.D0)

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST3=DIST**3

      CALL CALC_RR(RR,RW)

      A1=-P16*(U/DIST3-3.D0*RR/DIST3)

      RETURN
      END
      
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_RR_IJ(A1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST2,DIST3,RW(3),RR(3,3)
      REAL*8 aaI,aaJ
      REAL*8 M0RR
      PARAMETER(M0RR=1.D0/8.D0)
      REAL*8 PRE,PU,PRR,aaIaaJ,aaI2,aaJ2,aaI3,aaJ3

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_RR(RR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=M0RR/aaJ**3*U
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        aaI2=aaI**2
        aaJ2=aaJ**2
        aaI3=aaI**3
        aaJ3=aaJ**3
        PRE=1.D0/(8*64*aaI3*aaJ3)

        PU=  5*DIST3 - 27*DIST*(aaI2+aaJ2) 
     *     + 32*(aaI3+aaJ3)
     *     - 9*(aaI2-aaJ2)**2/DIST 
     *     - aaIaaJ**4*(aaI2 + 4*aaI*aaJ + aaJ2)/DIST3

        PRR=3*(aaIaaJ**2 - DIST2)**2
     *    *((aaI2 + 4*aaI*aaJ + aaJ2)-DIST2)/DIST3

        A1=PRE*(PU*U+PRR*RR)
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=M0RR/aaI**3*U
       ELSE
        DIST2=DIST**2
        DIST3=DIST**3
        aaIaaJ=aaI-aaJ
        aaI2=aaI**2
        aaJ2=aaJ**2
*-----------------------------------------------------------
        aaI3=aaI**3
        aaJ3=aaJ**3
*-----------------------------------------------------------
        PRE=1.D0/(8*64*(aaI**3)*(aaJ**3))

        PU=  5*DIST3 - 27*DIST*(aaI2+aaJ2) 
     *     + 32*(aaI3+aaJ3)
     *     - 9*(aaI2-aaJ2)**2/DIST 
     *     - aaIaaJ**4*(aaI2 + 4*aaI*aaJ + aaJ2)/DIST3

        PRR=3*(aaIaaJ**2 - DIST2)**2
     *    *((aaI2 + 4*aaI*aaJ + aaJ2)-DIST2)/DIST3

        A1=PRE*(PU*U+PRR*RR)
       ENDIF       
      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_RT_IJ(A1,R)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,DIST2,RW(3),EPSR(3,3)
      REAL*8 P18
      PARAMETER(P18=1.D0/8.D0)
      
      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST2=DIST**2
      
      CALL CALC_EPSR(EPSR,RW)

      A1=P18*EPSR/DIST2

      RETURN
      END
      
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_RT_IJ(A1,R,aaI,aaJ)
      IMPLICIT NONE
      REAL*8 A1(3,3),R(3)
      REAL*8 DIST,RW(3),EPSR(3,3)
      REAL*8 aaI,aaJ            ! radius
      REAL*8 PRE,PEPSR

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_EPSR(EPSR,RW)

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        A1=0.D0
       ELSE
        PRE=1.D0/(16*8*aaI**3*aaJ)
        PEPSR=((aaI-aaJ)+DIST)**2
     *   *(aaJ**2+2*aaJ*(aaI+DIST)-3*(aaI-DIST)**2)
     *   /DIST**2
        A1=PRE*PEPSR*EPSR
       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        A1=DIST/(8*(aaI**3))*EPSR
       ELSE
        PRE=1.D0/(16*8*aaI**3*aaJ)
        PEPSR=((aaI-aaJ)+DIST)**2
     *   *(aaJ**2+2*aaJ*(aaI+DIST)-3*(aaI-DIST)**2)
     *   /DIST**2
        A1=PRE*PEPSR*EPSR
       ENDIF       
      ENDIF
      
      RETURN
      END
      
C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_DD_IJ(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1(3,3,3,3),t4D1(3,3,3,3),t4D2(3,3,3,3),t4D0(3,3,3,3)
      REAL*8 DIST,DIST2,DIST5,RW(3),RR(3,3),R(3)
      REAL*8 aaI,aaJ,aaI2,aaJ2

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST5=DIST**5
      DIST2=DIST**2
      aaI2=aaI**2
      aaJ2=aaJ**2

      CALL CALC_RR(RR,RW)

      CALL CALC_t4D0(t4D0,RR)
      CALL CALC_t4D1(t4D1,RR)
      CALL CALC_t4D2(t4D2,RR)
      t4D2 = t4D2 - t4D1

      B1=0.D0

      B1=B1 +3.D0
     *      *(6.D0*(aaI2+aaJ2)-5.D0*DIST2)
     *      /(20.D0*DIST5)
     *      *t4D0

      B1=B1 -3.D0
     *      *(8.D0*(aaI2+aaJ2)-5.D0*DIST2)
     *      /(40.D0*DIST5)
     *      *t4D1

      B1=B1 +3.D0
     *      *(aaI2+aaJ2)
     *      /(20.D0*DIST5)
     *      *t4D2

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_DD_IJ_Y2(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1_CART(3,3,3,3),B1(5,5),R(3),a
      REAL*8 aaI,aaJ 
      INTEGER I,J

      CALL ROTNE_PRAGER_DD_IJ(B1_CART,R,aaI,aaJ)

      DO I=1,5
       DO J=1,5
        CALL mulT2aT4T2b(Y2(I,1:3,1:3),B1_CART,Y2(J,1:3,1:3),a) 
        B1(I,J)=a
       ENDDO
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_DD_IJ(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1(3,3,3,3),t4D1(3,3,3,3),t4D2(3,3,3,3),t4D0(3,3,3,3)
      REAL*8 DIST,DIST3,DIST5,RW(3),RR(3,3),R(3)
      REAL*8 aaI,aaJ,aaI2,aaJ2,aaI3,aaJ3

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST5=DIST**5
      DIST3=DIST**3
      aaI2=aaI**2
      aaJ2=aaJ**2
      aaI3=aaI**3
      aaJ3=aaJ**3

      CALL CALC_RR(RR,RW)

      CALL CALC_t4D0(t4D0,RR)
      CALL CALC_t4D1(t4D1,RR)
      CALL CALC_t4D2(t4D2,RR)
      t4D2 = t4D2 - t4D1

      B1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        B1= 3.D0/(aaJ**3*20.D0)*(t4D0+t4D1+t4D2)
       ELSE

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3)
     *     *(
     *     +3.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     -10.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3
     *     +32.D0*(aaI3+aaJ3)
     *     -30.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *      )*t4D0

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3)
     *     *(
     *     -2.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     +5.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3
     *     -15.D0*(aaI2-aaJ2)**2/DIST
     *     +32.D0*(aaI3+aaJ3)
     *     -25.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *     )*t4D1

      B1=B1+3.D0/(2560.D0*aaI3*aaJ3)
     *     *(
     *     +(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     -30.D0*(aaI2-aaJ2)**2/DIST
     *     +64.D0*(aaI3+aaJ3)
     *     -40.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *     )*t4D2

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        B1= 3.D0/(aaI**3*20.D0)*(t4D0+t4D1+t4D2)
       ELSE

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3)
     *     *(
     *     +3.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     -10.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3
     *     +32.D0*(aaI3+aaJ3)
     *     -30.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *      )*t4D0

      B1=B1+3.D0/(1280.D0*aaI3*aaJ3)
     *     *(
     *     -2.D0*(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     +5.D0*(aaI-aaJ)**4*(aaI2+4.D0*aaI*aaJ+aaJ2)/DIST3
     *     -15.D0*(aaI2-aaJ2)**2/DIST
     *     +32.D0*(aaI3+aaJ3)
     *     -25.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *     )*t4D1

      B1=B1+3.D0/(2560.D0*aaI3*aaJ3)
     *     *(
     *     +(aaI-aaJ)**6*(aaI2+6.D0*aaI*aaJ+aaJ2)/DIST5
     *     -30.D0*(aaI2-aaJ2)**2/DIST
     *     +64.D0*(aaI3+aaJ3)
     *     -40.D0*(aaI2+aaJ2)*DIST
     *     +5.D0*DIST3
     *     )*t4D2

       ENDIF       
      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_DD_IJ_Y2(B1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 B1_CART(3,3,3,3),B1(5,5),R(3),a
      REAL*8 aaI,aaJ
      INTEGER I,J

      CALL YAMAKAWA_DD_IJ(B1_CART,R,aaI,aaJ)

      DO I=1,5
       DO J=1,5
        CALL mulT2aT4T2b(Y2(I,1:3,1:3),B1_CART,Y2(J,1:3,1:3),a) 
        B1(I,J)=a
       ENDDO
      ENDDO

      RETURN
      END



C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_TD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3UR(3,3,3),t3RRR(3,3,3)
      REAL*8 DIST,DIST4,RW(3),R(3)
      REAL*8 aaI,aaJ,aaI2,aaJ2

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST4=DIST**4
      aaI2=aaI**2
      aaJ2=aaJ**2

      CALL CALC_UR(t3UR,RW)
      CALL CALC_RRR(t3RRR,RW)

      C1=0.D0

      C1= ( - 2.D0*(5.D0*aaI2*aaJ2+3.D0*aaJ2**2)/(5.D0*DIST4)*t3UR
     *      + aaJ2*(5.D0*aaI2+3.D0*aaJ2-3.D0*DIST**2)/DIST4*t3RRR
     *    )/(8.D0*aaJ2)

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_TD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ
      INTEGER I

      CALL ROTNE_PRAGER_TD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_TD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3UR(3,3,3),t3RRR(3,3,3)
      REAL*8 DIST,DIST2,DIST4,RW(3),R(3),PRE,PUR,PRRR
      REAL*8 aaI,aaJ

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST2=DIST**2
      DIST4=DIST**4

      CALL CALC_UR(t3UR,RW)
      CALL CALC_RRR(t3RRR,RW)

      C1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        C1=-3.D0*DIST/(20.D0*aaJ**3)*t3UR
       ELSE

        PRE = 1.D0/(8.D0*aaI*aaJ**3)

        PUR = (10.D0*DIST2 - 24.D0*aaI*DIST
     *       - 15.D0*(aaJ**2-aaI**2) 
     *       + (aaJ-aaI)**5*(aaI+5.D0*aaJ)/DIST4 )/40.D0

        PRRR = ((aaI-aaJ)**2-DIST2)**2
     *        *((aaI-aaJ)*(aaI+5*aaJ)-DIST2)
     *        /(16.D0*DIST4)

        C1 = PRE*(PUR*t3UR+PRRR*t3RRR)

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        C1=0.D0
       ELSE

        PRE = 1.D0/(8.D0*aaI*aaJ**3)

        PUR = (10.D0*DIST2 - 24.D0*aaI*DIST
     *       - 15.D0*(aaJ**2-aaI**2) 
     *       + (aaJ-aaI)**5*(aaI+5.D0*aaJ)/DIST4 )/40.D0

        PRRR = ((aaI-aaJ)**2-DIST2)**2
     *        *((aaI-aaJ)*(aaI+5*aaJ)-DIST2)
     *        /(16.D0*DIST4)

        C1 = PRE*(PUR*t3UR+PRRR*t3RRR)

       ENDIF       
      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_TD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ
      INTEGER I

      CALL YAMAKAWA_TD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_RD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3EPSRR(3,3,3)
      REAL*8 DIST,RW(3),R(3)
      REAL*8 aaI,aaJ

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST

      CALL CALC_EPSRR(t3EPSRR,RW)

      C1=0.D0

      C1= - 3.D0/(8.D0*DIST**3)*t3EPSRR

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE ROTNE_PRAGER_RD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ
      INTEGER I

      CALL ROTNE_PRAGER_RD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_RD_IJ(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1(3,3,3),t3EPSRR(3,3,3)
      REAL*8 DIST,DIST2,RW(3),R(3)
      REAL*8 aaI,aaJ

      DIST=SQRT( SUM(R**2) )
      RW=R/DIST
      DIST2=DIST**2

      CALL CALC_EPSRR(t3EPSRR,RW)

      C1=0.D0

      IF (aaI.LE.aaJ) THEN
       IF (DIST.LE.(aaJ-aaI)) THEN
        C1=0.D0
       ELSE

        C1 = -15.D0/(1280.D0*(aaI**3)*(aaJ**3)*DIST**3)
     *       *((aaI-aaJ)**2 - DIST2)**2
     *       *((aaI**2 + 4*aaI*aaJ + aaJ**2)-DIST2)
     *       *t3EPSRR

       ENDIF
      ELSE
       IF (DIST.LE.(aaI-aaJ)) THEN
        C1=0.D0
       ELSE

        C1 = -15.D0/(1280.D0*(aaI**3)*(aaJ**3)*DIST**3)
     *       *((aaI-aaJ)**2 - DIST2)**2
     *       *((aaI**2 + 4*aaI*aaJ + aaJ**2)-DIST2)
     *       *t3EPSRR

       ENDIF       
      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE YAMAKAWA_RD_IJ_Y2(C1,R,aaI,aaJ)
      USE TENSORS
      IMPLICIT NONE
      REAL*8 C1_CART(3,3,3),C1(3,5),R(3),V(3)
      REAL*8 aaI,aaJ
      INTEGER I

      CALL YAMAKAWA_RD_IJ(C1_CART,R,aaI,aaJ)

      DO I=1,5
       CALL mulT3T2(C1_CART,Y2(I,1:3,1:3),V)
       C1(1:3,I)=V
      ENDDO

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE MOB_TO_FRI(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN)

      CALL MATREV(APP,6*NN)

      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ + MATMUL(AQP,APQ)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE FRI_TO_MOB(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN)

      CALL MATREV(APP,6*NN)

      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ - MATMUL(AQP,APQ)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE FRI_TO_MOB_RED(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN)

      AQP = TRANSPOSE(APQ)
      APQ = MATMUL(APP,APQ)
      AQQ = AQQ - MATMUL(AQP,APQ)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE INVFRI_TO_FRI(APP,APQ,AQQ,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 APP(6*NN,6*NN),APQ(6*NN,5*NN),AQQ(5*NN,5*NN)
      REAL*8 AQP(5*NN,6*NN),AP(11*NN,11*NN)

      AP(11*NN,11*NN) = 0.D0

      AP(1:6*NN,1:6*NN) = APP
      AP(1:6*NN,6*NN+1:11*NN) = APQ 
      AP(6*NN+1:11*NN,1:6*NN) = TRANSPOSE(APQ)
      AP(6*NN+1:11*NN,6*NN+1:11*NN) = AQQ 

      CALL MATREV(APP,6*NN)
      CALL MATREV(AP,11*NN)

      APP = AP(1:6*NN,1:6*NN)
      APQ = AP(1:6*NN,6*NN+1:11*NN)
      AQQ = AP(6*NN+1:11*NN,6*NN+1:11*NN)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE MATREV(A,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 A(NN,NN)
      INTEGER INFO,I,J
      CHARACTER*1 UPLO
      PARAMETER (UPLO = 'U')

      CALL DPOTRF(UPLO,NN,A,NN,INFO)
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
      IF(INFO.GT.0) THEN
      WRITE(*,*) 'k=',INFO
      WRITE(*,*) 'The leading minor of order k is not'
      WRITE(*,*) 'positive definite, and the factorization could not be'
      WRITE(*,*) 'completed.'
      ENDIF
      CALL DPOTRI(UPLO,NN,A,NN,INFO)

      DO I=2,NN
       DO J=1,I-1
        A(I,J)=A(J,I)
       ENDDO
      ENDDO

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************


C May 24, 2003
C CHOLESKY FROM LAPACK.

      SUBROUTINE CHOLESKY(B,NN)
      IMPLICIT NONE
      INTEGER NN
      REAL*8 B(NN,NN)
      CHARACTER*1 UPLO
      PARAMETER (UPLO = 'L')
      INTEGER I,J,INFO

      CALL DPOTRF(UPLO,NN,B,NN,INFO)
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
      IF(INFO.GT.0) THEN
      WRITE(*,*) 'k=',INFO
      WRITE(*,*) 'The leading minor of order k is not'
      WRITE(*,*) 'positive definite, and the factorization could not be'
      WRITE(*,*) 'completed.'
      STOP
      ENDIF

C     RETURN     !!!  The strictly upper triangular part of A is not referenced.
      DO I=1,NN-1
      DO J=I+1,NN
       B(I,J)=0.D0
      ENDDO
      ENDDO

      RETURN
      END

C***************************************************************************
C***************************************************************************
C***************************************************************************

