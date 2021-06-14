C WRITES THE CURRENT CONFIGURATION ON FILE

      PROGRAM MAIN
      USE SIZE           ! NN, MM, NSpP 
      USE CONFIG         ! CONF,RADII,POLY_LEN
      USE FORCE_PAR      ! K_PAR, ELONGMAX 
      USE LATTICE_SKEW   ! ERR
      USE STEPPER_PAR    ! DT, SEED
      USE TENSORS    ! DT, SEED
      USE MULT       ! LL
      IMPLICIT NONE
      INTEGER I,J, SEEDSIZE, INP, K,K_WRITE
      INTEGER,ALLOCATABLE :: SEED(:)
      LOGICAL THERE
      REAL*8 T, DR2, DR(3),DRO(3), END_TIME, BASE_DT,
     +       EXPER,TMPR(3),DROLD(3)
      REAL*8 PHI,VOL


************ PARAMS *******************************

      BASE_DT  = 0.1D0
      END_TIME = 10000.D0
      K_WRITE = 1
      DRO = 0.D0

*--------------------------------

      CALL INIT_U
      CALL INIT_EPS
      CALL INIT_Y2

C     program arguments
C     NN, END_TIME, BASE_DT, K_WRITE, K_PAR, ELONGMAX
C     which HAVE to be specified

      INQUIRE( FILE='control_file.dat', EXIST=THERE ) 
       IF ( THERE ) THEN
       OPEN(10,FILE='control_file.dat')
       READ(10,*) NN
       READ(10,*) END_TIME 
       READ(10,*) BASE_DT
       READ(10,*) K_WRITE
       DT = BASE_DT
       CALL RANDOM_SEED(SIZE=SEEDSIZE)
       ALLOCATE(SEED(SEEDSIZE))
       READ(10,*) INP
       SEED = INP
       CALL RANDOM_SEED(PUT=SEED)
       READ(10,*) K_PAR
       READ(10,*) ELONGMAX
       CLOSE(10)
      ELSE
       WRITE(0,*) "error: no control_file.dat file "
     *   // "present in this folder"
       WRITE(0,*) "        Please provide one"
       CALL EXIT()
      ENDIF
      
      WRITE(*,*) 'N: ',NN,' dt: ',BASE_DT
      DROLD = 0.D0
      K = 0

      CALL CONFIGURATION
      CALL POS_CORRECTION(CONF,RADII,DRO)

      CALL INIT_VMD_WRITE(RADII)
      CALL INIT_CM_WRITE()

      DO
       K=K+1

       CALL POS_CORRECTION(CONF,RADII,DR)

       IF (mod(K,K_WRITE) .eq. 0) then
        CALL VMD_WRITE(CONF,T,K)
        DR2 = SQRT(SUM((DR-DRO)*(DR-DRO)))
        CALL CM_WRITE(K,T,DR2,DR)
        CALL PROGRESSSTATUS(INT(T/END_TIME*100))
       ENDIF

       CALL STEPPER(CONF,RADII,DT,T)

       IF (T.GT.END_TIME) THEN
        CALL PROGRESSSTATUS(INT(100))
        WRITE(*,*)
        EXIT
       ENDIF
          
      ENDDO

      CALL END_VMD_WRITE()
      CALL END_CM_WRITE()

      STOP
      END


