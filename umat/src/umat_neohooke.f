C=======================================================================
C     UMAT: compressible Neo-Hookean hyperelasticity, Abaqus form
C
C       U = C10 (I1bar - 3) + (J - 1)**2 / D1
C
C     PROPS(1) = C10 (stress), PROPS(2) = D1 (1/stress)
C     Small-strain limit: shear modulus 2*C10, bulk modulus 2/D1.
C     Matches Abaqus/CalculiX *HYPERELASTIC, NEO HOOKE with the same
C     C10 and D1. Requires NLGEOM. No state variables.
C
C     Output frame
C       Abaqus (default): Cauchy stress in the global basis, built
C         from Bbar = J**(-2/3) F F^T.
C       CalculiX (compile with -DCCX_COROT, material name ABAQUSNL...):
C         ccx's umat_abaqusnl wrapper expects the corotated Cauchy
C         stress R^T sigma R, which for an isotropic material is the
C         same formula with Cbar = J**(-2/3) F^T F in place of Bbar.
C
C     DDSDDE is (1/J) times the Jaumann tangent of the Kirchhoff
C     stress, the tangent Abaqus/Standard expects for finite-strain
C     UMATs; engineering shear strains in columns 4-6.
C
C     Voigt order 11, 22, 33, 12, 13, 23.
C=======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(*),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C
      DIMENSION BB(3,3), DLT(3,3), TAU(3,3)
      INTEGER IV(6), JV(6)
      DATA IV /1,2,3,1,1,2/, JV /1,2,3,2,3,3/
C
      C10 = PROPS(1)
      D1  = PROPS(2)
C
C     Jacobian J = det F
      DET = DFGRD1(1,1)*(DFGRD1(2,2)*DFGRD1(3,3)
     1                  -DFGRD1(2,3)*DFGRD1(3,2))
     2    - DFGRD1(1,2)*(DFGRD1(2,1)*DFGRD1(3,3)
     3                  -DFGRD1(2,3)*DFGRD1(3,1))
     4    + DFGRD1(1,3)*(DFGRD1(2,1)*DFGRD1(3,2)
     5                  -DFGRD1(2,2)*DFGRD1(3,1))
      IF (DET .LE. 0.D0) THEN
        PNEWDT = 0.25D0
        RETURN
      END IF
      SCALE = DET**(-2.D0/3.D0)
C
C     Isochoric left (Abaqus) or right (CalculiX) Cauchy-Green tensor
      DO I = 1, 3
        DO J = 1, 3
          SUM = 0.D0
          DO K = 1, 3
#ifdef CCX_COROT
            SUM = SUM + DFGRD1(K,I)*DFGRD1(K,J)
#else
            SUM = SUM + DFGRD1(I,K)*DFGRD1(J,K)
#endif
          END DO
          BB(I,J) = SCALE*SUM
          DLT(I,J) = 0.D0
        END DO
        DLT(I,I) = 1.D0
      END DO
      TRB = BB(1,1) + BB(2,2) + BB(3,3)
C
C     Kirchhoff stress, then Cauchy = Kirchhoff / J
      PVOL = 2.D0/D1*DET*(DET - 1.D0)
      DO I = 1, 3
        DO J = 1, 3
          TAU(I,J) = 2.D0*C10*(BB(I,J) - TRB/3.D0*DLT(I,J))
     1             + PVOL*DLT(I,J)
        END DO
      END DO
      DO M = 1, NTENS
        STRESS(M) = TAU(IV(M),JV(M))/DET
      END DO
C
C     Tangent: (1/J) d(tau)/d(D), Jaumann
      EVOL = 2.D0/D1*DET*(2.D0*DET - 1.D0)
      DO M = 1, NTENS
        I = IV(M)
        J = JV(M)
        DO N = 1, NTENS
          K = IV(N)
          L = JV(N)
          DDSDDE(M,N) = (2.D0*C10*(
     1        0.5D0*(DLT(I,K)*BB(J,L) + BB(I,K)*DLT(J,L)
     2             + DLT(I,L)*BB(J,K) + BB(I,L)*DLT(J,K))
     3      - 2.D0/3.D0*(BB(I,J)*DLT(K,L) + DLT(I,J)*BB(K,L))
     4      + 2.D0/9.D0*TRB*DLT(I,J)*DLT(K,L))
     5      + EVOL*DLT(I,J)*DLT(K,L))/DET
        END DO
      END DO
C
C     Elastic strain energy density
      SSE = C10*(TRB - 3.D0) + (DET - 1.D0)**2/D1
      RETURN
      END
