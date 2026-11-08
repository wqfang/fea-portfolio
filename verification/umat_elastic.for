C     Minimal isotropic linear-elastic UMAT, used only to test that
C     Abaqus can compile and link user subroutines.
C     PROPS(1) = Young's modulus E, PROPS(2) = Poisson's ratio nu
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C
      E     = PROPS(1)
      ANU   = PROPS(2)
      ALAMB = E*ANU/((1.D0+ANU)*(1.D0-2.D0*ANU))
      AMU   = E/(2.D0*(1.D0+ANU))
C
C     Elastic stiffness (Abaqus uses engineering shear strain)
      DO I = 1, NTENS
        DO J = 1, NTENS
          DDSDDE(I,J) = 0.D0
        END DO
      END DO
      DO I = 1, NDI
        DO J = 1, NDI
          DDSDDE(I,J) = ALAMB
        END DO
        DDSDDE(I,I) = ALAMB + 2.D0*AMU
      END DO
      DO I = NDI+1, NTENS
        DDSDDE(I,I) = AMU
      END DO
C
C     Stress update: sigma = sigma + C : d(epsilon)
      DO I = 1, NTENS
        DO J = 1, NTENS
          STRESS(I) = STRESS(I) + DDSDDE(I,J)*DSTRAN(J)
        END DO
      END DO
C
      RETURN
      END
