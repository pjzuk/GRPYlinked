      SUBROUTINE FORCE(FEXT,CONFL)
      USE SIZE        ! NN
      USE CONFIG      ! 
      USE FORCE_PAR   ! K_PAR, ELONGMAX
      IMPLICIT NONE
      REAL*8 FEXT(3*NN),CONFL(3,NN)
      REAL*8 FS(3*NN)   ! STRECHING FORCES
      REAL*8 T(3,0:NN+2)
      REAL*8 L(1:NN+1),LINV(1:NN+1)
      REAL*8 R,L0,T0(3)
      INTEGER I,J,K,M
      REAL*8 B1(3),B2(3)
      REAL*8 CS,S,TH,DF
      REAL*8 a11,a12,a22

      ELONGMAX = 0.1D0

      T=0.D0
      L   =0.D0
      LINV=0.D0
      DO I=1,NN-1
       T(:,I)=CONFL(:,I)-CONFL(:,I+1)
       L(I)=SQRT(SUM(T(:,I)*T(:,I)))      
        IF (L(I).GT.0.D0) THEN
         T(:,I)=T(:,I)/L(I)
        ENDIF
         LINV(I)=1.D0/L(I)
      ENDDO

C STRETCHING FORCE FENE -----------------------------------------

      FS=0.D0
      DO K=1,NN-1
       IF (L(K).GE.LINK(K)) THEN
        FS(3*(K-1)+1:3*(K-1)+3) = FS(3*(K-1)+1:3*(K-1)+3)
     * -(L(K)-LINK(K))/(1-((L(K)/LINK(K)-1)/ELONGMAX)**2)*T(:,K)
        FS(3*K+1:3*K+3) = FS(3*K+1:3*K+3) 
     * +(L(K)-LINK(K))/(1-((L(K)/LINK(K)-1)/ELONGMAX)**2)*T(:,K)
       ENDIF
       IF ((L(K)-LINK(K)).GT.(.5D0)) THEN
        WRITE (*,*) "LINK ",K-1," TO ",K," GOT DISCONECTED"
        STOP
       ENDIF
      ENDDO
      FS=K_PAR*FS

CC TOTAL FORCE -------------------------------------------------

      FEXT = FS

CC TOTAL FORCE SET TO 0 ----------------------------------------
C
C      FEXT=0.D0

      RETURN
      END

