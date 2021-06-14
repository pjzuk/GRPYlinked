      SUBROUTINE CONFIGURATION
      USE SIZE           ! NN, MM, NSpP
      USE CONFIG         ! CONF,RADII
      USE FORCE_PAR
      IMPLICIT NONE
      REAL*8 X,Y,Z,L,R(3),PI,L1_PAR
      INTEGER I,J,K,M
      LOGICAL THERE
      PARAMETER(PI=3.141592653589793D0)
      
      IF(.NOT.ALLOCATED(CONF)) ALLOCATE( CONF(3,NN) )
      IF(.NOT.ALLOCATED(RADII)) ALLOCATE( RADII(NN) )
      IF(.NOT.ALLOCATED(LINK)) ALLOCATE( LINK(NN-1) )

      L1_PAR = NN**(1.D0/3.D0)*2.D0/3.D0

      INQUIRE( FILE='model.dat', EXIST=THERE ) 
      IF ( THERE ) THEN
       OPEN(11,FILE='model.dat')
       DO I=1,NN
        READ(11,*) RADII(I)
       ENDDO
       DO I=1,NN-1
        READ(11,*) LINK(I)
       ENDDO
       CLOSE(11)
      ELSE
       WRITE(0,*) "error: no model.dat file "
     *   // "present in this folder"
       WRITE(0,*) "        Please provide one"
       CALL EXIT()
      ENDIF

      CONF(1,1) = 0.D0 
      CONF(2,1) = 0.D0 
      CONF(3,1) = 0.D0 

      X=1.D0
      Y=0.D0
      Z=0.D0

      DO I = 2,NN

       CONF(1,I) = CONF(1,(I-1)) +
     *        LINK(I-1)*X 
       CONF(2,I) = CONF(2,(I-1)) +
     *        LINK(I-1)*Y
       CONF(3,I) = CONF(3,(I-1)) +
     *        LINK(I-1)*Z

      ENDDO    

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE POS_CORRECTION(CONF,RADII,DR)
      USE SIZE
      IMPLICIT NONE
      REAL*8 DR(3),CONF(3,NN),RADII(NN),MASS
      INTEGER I

C   CORRECTIONS FOR CM DISPLACEMENT

      MASS = SUM(RADII**3.D0)
      DR = 0.D0
      DO I = 1,NN
       DR = DR + CONF(:,I)*(RADII(I)**3.D0)
      ENDDO
      DR = DR/MASS

      RETURN
      END
