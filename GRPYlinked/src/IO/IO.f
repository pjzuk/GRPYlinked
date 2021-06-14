
      SUBROUTINE PROGRESSSTATUS(PROC)
      IMPLICIT NONE
      INTEGER PROC
      CHARACTER*10 :: BAR = '  0% DONE '

      WRITE(UNIT=BAR(1:3),FMT='(i3)') PROC

      WRITE(*,FMT='(A,A10,$)') CHAR(13),BAR

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

C VTF file for vmd reader header

      SUBROUTINE INIT_VMD_WRITE(RADII)
      USE SIZE
      USE FORCE_PAR
      USE STEPPER_PAR
      IMPLICIT NONE
      REAL(8)  RADII(NN)
      INTEGER I

      OPEN(31,FILE='out_trajectory.VTF') ! xyz file for vmd

      WRITE(31,*) '# ',DT
      DO I=1,NN
       WRITE(31,*) 'atom ',I-1,'radius',RADII(I),'name C'
      ENDDO

      WRITE(31,*) ''

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE VMD_WRITE(CONF,T,K)
      USE SIZE
      IMPLICIT NONE
      REAL(8) CONF(3,NN),CONFQ(0:3,NN),T,DR(3),ROTM(3,3)
      INTEGER I,K

      WRITE(31,*) ""
      WRITE(31,*) "timestep"
      WRITE(31,*) "#",T,K
      DO I = 1, NN
       WRITE(31,*) CONF(:,I)
      ENDDO

      CALL FLUSH(31)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

C close VTF file for vmd reader

      SUBROUTINE END_VMD_WRITE()
      IMPLICIT NONE

      CLOSE(31) ! xyz file for vmd

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

      SUBROUTINE INIT_CM_WRITE()
      IMPLICIT NONE

      OPEN(32,FILE='out_cm.dat') 

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************


      SUBROUTINE CM_WRITE(K,T,DR2,DR)
      USE SIZE
      IMPLICIT NONE
      REAL(8) T,DR2,DR(3)
      INTEGER K

      WRITE (32,*) K,T,DR2,DR
      CALL FLUSH(32)

      RETURN
      END


C***********************************************************
C***********************************************************
C***********************************************************

C close VTF file for vmd reader

      SUBROUTINE END_CM_WRITE()
      IMPLICIT NONE

      CLOSE(32) 

      RETURN
      END

C***********************************************************
C***********************************************************
C***********************************************************

