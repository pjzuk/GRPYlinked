*****************************************

      MODULE SIZE       ! NN,NN6,LLG
      INTEGER NN,NN6,LLG
      END MODULE SIZE

*****************************************

      MODULE MULT       ! LL
      INTEGER :: LL=2
      END MODULE MULT

*****************************************

      MODULE LATTICE_SKEW 
      REAL*8 LR(3,3),LI(3,3),LV
      REAL*8 SIG,ERR,LOGERR
      INTEGER NR,NI
      END MODULE LATTICE_SKEW

*****************************************

      MODULE MULT_FULL  ! LL2,LLF,LLM,LLG1
      INTEGER LL2,LLF,LLM,LLG1
      END MODULE MULT_FULL

*****************************************

      MODULE Z_COEF_polyd
      REAL*8, ALLOCATABLE :: Z0(:,:),Z1(:,:),Z2(:,:),ZB(:,:)
      END MODULE Z_COEF_polyd

*****************************************

      MODULE S_COEF
      REAL*8, ALLOCATABLE :: S11P(:,:),So00P(:,:),S00P(:,:),S01P(:,:),
     *                       S02P(:,:),S10P(:,:),S20P(:,:)
      REAL*8, ALLOCATABLE :: S11N(:,:),So00N(:,:),S00N(:,:),S01N(:,:),
     *                       S02N(:,:),S10N(:,:),S20N(:,:)
      END MODULE S_COEF

*****************************************

      MODULE FACT
      REAL*8, ALLOCATABLE :: SIL(:)
      END MODULE FACT

*****************************************

      MODULE GAM0_MOD     !  GAM0
      REAL*8, ALLOCATABLE :: GAM0(:)
      END MODULE GAM0_MOD

*****************************************

      MODULE PI_I_MOD   !  PI_I1,PI_I2
      COMPLEX*16, ALLOCATABLE :: PI_I1(:),PI_I2(:)
      END MODULE PI_I_MOD

*****************************************

      MODULE THREADS
      INTEGER NCORES
      END MODULE THREADS

*****************************************

