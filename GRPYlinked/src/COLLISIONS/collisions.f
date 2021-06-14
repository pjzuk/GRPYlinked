C COLLISIONS MANAGEMENT ----------------------------------------
C THIS PROCEDURES TAKES AN INITIAL CONFIGURATION, DISPLACEMENMT, RADII,
C DT AND T AND EXECUTES COLLISIONS THAT OCCUR FROM T TO T+DT
C AFTER EACH COLLISION THE LATTICE SKEW IS CORRECTED ACCORDING TO ELAPSED TIME
C IT IS REPEATED UNTIL TIME EXEEDS DT OR NO MORE COLLISIONS WITHIN REMAINING TIME

      SUBROUTINE COLLISIONS_PER(CONF,DISP,RADII,DT_ORG,T_ORG)
      USE SIZE        ! NN
      IMPLICIT NONE
      REAL*8 CONF(3*NN),DISP(3*NN),RADII(NN),DT_ORG,T_ORG,T
      REAL*8 VEL(3*NN),RADIITWO(2),VEL_2(6),CONF_2(6)
      REAL*8 DTCOL,DT
      INTEGER I,J,COL_N,kk
      
      DT=DT_ORG
      T=T_ORG
      DO I=1,3*NN
        IF (DISP(I).EQ.0) THEN
          VEL=0
        ELSE
          VEL=DISP/DT
        ENDIF 
      ENDDO

      kk=0  ! number of collisions
      DO

      CALL CHECK_ALL_FOR_COLL(COL_N,DTCOL,I,J,CONF,VEL,RADII)

      IF(COL_N==0) THEN
       CONF=CONF+VEL*DT
       EXIT
      ENDIF

      IF (DT<DTCOL) THEN
       CONF=CONF+VEL*DT
       EXIT
      ENDIF

      kk=kk+1
      RADIITWO(1)=RADII(I)
      RADIITWO(2)=RADII(J)

C THIS IS EXECUTION OF COLLISION ------- 
      VEL_2(1:3)=VEL(3*(I-1)+1:3*I)
      VEL_2(4:6)=VEL(3*(J-1)+1:3*J)
      CONF_2(1:3)=CONF(3*(I-1)+1:3*I)
      CONF_2(4:6)=CONF(3*(J-1)+1:3*J)

      CALL COLLISION_EXEC(CONF_2,VEL_2,RADIITWO,DTCOL)
      CONF=CONF+VEL*DTCOL
      CONF(3*(I-1)+1:3*I)=CONF_2(1:3)
      CONF(3*(J-1)+1:3*J)=CONF_2(4:6)
      VEL(3*(I-1)+1:3*I)=VEL_2(1:3)
      VEL(3*(J-1)+1:3*J)=VEL_2(4:6)
      DT=DT-DTCOL
      T = T + DTCOL
        
      ENDDO

C      write(22,*) kk,(1.D0*kk)/NN
C      backspace(22)
C      read(22,*)
  
      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

C CHECK ALL PARTICLES FOR COLLISONS ---------------------------
C THIS PROCEDURE CHECKS ALL PARTICLES IN NAIVE WAY FOR COLLISIONS

      SUBROUTINE CHECK_ALL_FOR_COLL(COL_N,DTCOL,I,J,CONF,VEL,RADII)
      USE SIZE
      USE FORCE_PAR
      IMPLICIT NONE
      INTEGER COL_N,I,J
      REAL*8 DTCOL,CONF(3*NN),VEL(3*NN),RADII(NN)
      REAL*8 CONF_REL(3),VEL_REL(3),RADIITWO(2),DIST,CONF_REL_Z
      REAL*8 CONFJO(3),VELJO(3)
      REAL*8 COLL_TIME(NN*(NN+1)/2)
      INTEGER COLL_IJ(2,NN*(NN+1)/2),K
      LOGICAL IF_COL

      COL_N=0

      DO I=1,NN-1
       DO J=I+1,NN
        CONFJO=CONF(3*(J-1)+1:3*J)
        VELJO = VEL(3*(J-1)+1:3*J)
        CONF_REL   = CONF(3*(J-1)+1:3*J) - CONF(3*(I-1)+1:3*I)
        CONF_REL_Z = CONF_REL(3)
        DIST=SQRT(SUM(CONF_REL*CONF_REL))
        CONF(3*(J-1)+1:3*J) = CONF(3*(I-1)+1:3*I) + CONF_REL

        VEL(3*(J-1)+1)=VEL(3*(J-1)+1)-(CONF_REL_Z - CONF_REL(3))*GAMMA

        VEL_REL=VEL(3*(J-1)+1:3*J) - VEL(3*(I-1)+1:3*I)

        RADIITWO(1) = RADII(I)
        RADIITWO(2) = RADII(J)
        CALL COLLISION_CHECK(IF_COL,CONF_REL,VEL_REL,RADIITWO,DTCOL) 
        IF (IF_COL) THEN
         COL_N = COL_N + 1
         COLL_TIME(COL_N) = DTCOL
         COLL_IJ(1,COL_N) = I
         COLL_IJ(2,COL_N) = J
        ENDIF
        CONF(3*(J-1)+1:3*J)=CONFJO
        VEL(3*(J-1)+1:3*J)=VELJO
       ENDDO
      ENDDO

      IF(COL_N/=0) THEN
       K=MINLOC(COLL_TIME(1:COL_N),1)
       I=COLL_IJ(1,K)
       J=COLL_IJ(2,K)
       DTCOL=COLL_TIME(K)
      ENDIF

      RETURN
      END
 
C***********************************************************
C***********************************************************
C***********************************************************

C COLLISION DETECTION -----------------------------------------
C IF PARTICLES WITH GIVEN POSITIONS AND VELOCITY COLLIDE WITHIN DT
C PROCEDURE RETURNS LOGICAL IF_COL TRUE AND DT AS COLLISION TIME

      SUBROUTINE COLLISION_CHECK(IF_COL,CONF_REL,VEL_REL,RADIITWO,DT)
      IMPLICIT NONE
      REAL*8 CONF_REL(3),VEL_REL(3),RADIITWO(2)
      REAL*8 DT,DELTA,TIME,B,SQRTDELTA
      LOGICAL IF_COL

      B = DOT_PRODUCT(CONF_REL,VEL_REL)

      IF (B.GT.0.D0) THEN
       IF_COL = .FALSE.
       RETURN
      ENDIF

      DELTA = B**2
     * - DOT_PRODUCT(VEL_REL,VEL_REL)*(DOT_PRODUCT(CONF_REL,CONF_REL)
     *    - (RADIITWO(1) + RADIITWO(2))**2)

      IF (DELTA < 0.D0) THEN
       IF_COL = .FALSE.
      ELSE
       SQRTDELTA = SQRT(DELTA)
       TIME = MIN(MAX(-B - SQRTDELTA,0.D0),MAX(-B + SQRTDELTA,0.D0))
       IF (TIME.GT.0D0) THEN
        TIME = TIME/DOT_PRODUCT(VEL_REL,VEL_REL)
        DT=TIME
        IF_COL = .TRUE.
       ELSE
        IF_COL = .FALSE.
       ENDIF
      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

C COLLISION EXECUTION -----------------------------------------
C GIVEN THE COLLISION TIME PROCEDURE SETS PARTICLE IN CONTACT WITH APPROPRIATE VELOCITIES
C THAT ARE OBTAINED JUST AFTER COLLISION
C MASS OF THE PARTICLES IS CALCULATED AS RADII**3

      SUBROUTINE COLLISION_EXEC(CONF_2,VEL_2,RADIITWO,DT)
      USE FORCE_PAR
      IMPLICIT NONE
      REAL*8 CONF_2(6),VEL_2(6),RADIITWO(2),CONF_COL(3),VEL_CORR
      REAL*8 DT,M1,M2,VEL_RAD_2(2),CONF_COL_O(3),DIST
      REAL*8 VEL_CM(3)

      CONF_COL_O = CONF_2(4:6) - CONF_2(1:3)
      CONF_COL=CONF_COL_O
      DIST = SQRT(SUM(CONF_COL*CONF_COL))

      M1 = RADIITWO(1)**3
      M2 = RADIITWO(2)**3

      IF ((ABS(SUM(CONF_COL_O**2) - DIST**2)/SUM(CONF_COL_O**2))
     * .LE.1D-10) THEN
       CONF_2 = CONF_2 + VEL_2*DT
       VEL_CM = (M1*VEL_2(1:3) + M2*VEL_2(4:6))/(M1+M2)
       VEL_2(1:3) = VEL_2(1:3) - VEL_CM
       VEL_2(4:6) = VEL_2(4:6) - VEL_CM
       CONF_COL = CONF_2(4:6) - CONF_2(1:3)
       CONF_COL = CONF_COL/SQRT(SUM(CONF_COL**2))
       VEL_RAD_2(1)=DOT_PRODUCT(VEL_2(1:3),CONF_COL)
       VEL_RAD_2(2)=DOT_PRODUCT(VEL_2(4:6),CONF_COL)

       VEL_2(1:3) = VEL_2(1:3)
     *  + 2.D0*M2*(VEL_RAD_2(2) - VEL_RAD_2(1))
     *           /(M2+M1)*CONF_COL
       VEL_2(4:6) = VEL_2(4:6)
     *  + 2.D0*M1*(VEL_RAD_2(1) - VEL_RAD_2(2))
     *           /(M2+M1)*CONF_COL

       VEL_2(1:3) = VEL_2(1:3) + VEL_CM
       VEL_2(4:6) = VEL_2(4:6) + VEL_CM

      ELSE

       VEL_CORR=(CONF_COL(3)-CONF_COL_O(3))*GAMMA
       CONF_2(4:6) = CONF_2(1:3) + CONF_COL
       VEL_2(4) = VEL_2(4) + VEL_CORR

       CONF_2(1:3) = CONF_2(1:3) + VEL_2(1:3)*DT
       CONF_2(4:6) = CONF_2(4:6) + VEL_2(4:6)*DT

       VEL_CM = (M1*VEL_2(1:3) + M2*VEL_2(4:6))/(M1+M2)
       VEL_2(1:3) = VEL_2(1:3) - VEL_CM
       VEL_2(4:6) = VEL_2(4:6) - VEL_CM

       CONF_COL = CONF_2(4:6) - CONF_2(1:3)
       CONF_COL = CONF_COL/SQRT(SUM(CONF_COL**2))

       VEL_RAD_2(1)=DOT_PRODUCT(VEL_2(1:3),CONF_COL)
       VEL_RAD_2(2)=DOT_PRODUCT(VEL_2(4:6),CONF_COL)

       VEL_2(1:3) = VEL_2(1:3)
     *  + 2.D0*M2*(VEL_RAD_2(2) - VEL_RAD_2(1))
     *           /(M2+M1)*CONF_COL
       VEL_2(4:6) = VEL_2(4:6)
     *  + 2.D0*M1*(VEL_RAD_2(1) - VEL_RAD_2(2))
     *           /(M2+M1)*CONF_COL

       VEL_2(1:3) = VEL_2(1:3) + VEL_CM
       VEL_2(4:6) = VEL_2(4:6) + VEL_CM

      ENDIF

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************


