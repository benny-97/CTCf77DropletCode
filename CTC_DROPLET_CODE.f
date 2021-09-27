C*********************************************************************
C                                                                    *
C      PROGRAM FOR QUASI-ONE-DIMENSIONAL COMPRESSIBLE FLOW IN A      *
C      VARIABLE DUCT USING CONSERVATIVE VARIABLES. IT IS ASSUMED     *
C      THAT THE GAS IS IDEAL AND THAT THE VISCOSITY VARIES WITH THE  *
C      SQUARE ROOT OF THE TEMPERATURE.                               *
C                                                                    *
C*********************************************************************
C2345678
      DIMENSION XL(101),DIA(101),AR(101),U(101),P(101),T(101),IC(101)
      DIMENSION DEN(101),UP(101),TP(101)
      COMMON TAUF,VISC,AG,FRICP,DX,DIAAV,PRN,SHR,UAVE,DENAV,DP,TMELT,
     1       TW,RLC,RADF,X1,CPS,TAVE,ZL
C
C---------------------------------------------------------------------
C
C      INPUT PARAMETERS
C---------------------------------------------------------------------
       OPEN(UNIT=1,FILE='nozzpar.dat',STATUS='UNKNOWN')
C--INLET DIAMETER (M)
       DIA1=0.05
C--EXIT DIAMETER (M)
       DIA2=0.025
C--NOZZLE LENGTH (M)
       AL=0.5
C--PIPE ROUGHNESS (M)
       ROUGH=0.003
C--INLET PRESSURE (PA)
       P0=101.E03
C--INLET TEMPERATURE (K)
       T0=293.0
C--INLET VELOCITY (M/S)
       U0=18.0
C--GAS CONSTANT (J/KG-K)
       RG=287
C--RATIO OF SPECIFIC HEATS
       AK=1.4
C--VISCOSITY (N-2/M2)
       VISCO=3.0E-05
C--PRANDTL NUMBER
       PRN=0.72
C--WALL TEMPERATURE (K)
       TW=293.0
C---------------------------------------------------------------------
C
C      PARTICLE PARAMETERS
C---------------------------------------------------------------------
C--PARTICLE DIAMETER (meters)
       DP=100.0E-06
C--PARTICLE DENSITY (kg/m3)
       RHOP=2500.0
C--LOADING
       ZL=1.0
C--PARTICLE WALL FRICTION COEFFICIENT
       FRICP=0.001
C--ACCELERATION DUE TO GRAVITY (m/s2)
       AG=0.0
C--INITIAL PARTICLE VELOCITY (m/s)
       UP0=18.0
C--INITIAL TEMPERATURE
       TP0=1000.0
C--PARTICLE SPECIFIC HEAT (J/kg/K)
       CPS=800.0
C--LATENT HEAT OF FUSION (J/kg)
       HLAT=10000.0
C--MELTING TEMPERATURE (K)
       TMELT=1200.0
C--EMISSIVITY OF PARTICLE SURFACE
       EPSILON=0.0
C---------------------------------------------------------------------
C
C      COMPUTATIONAL PARAMETERS
C---------------------------------------------------------------------
C--NUMBER OF GRID LINES
       IMAX=101
C--LINES SKIPPED IN PRINTOUT
       ISKIP=2
C--RESIDUAL FOR VELOCITY CONVERGENCE
       RESMAX=0.001
C---------------------------------------------------------------------
       PLANK=5.668E-08
       RES=RESMAX*U0
       XL(1)=0.0
       DX=AL/(IMAX-1)
       FAC1=AK/(AK+1)
       FAC2=2*(AK*AK-1)/AK/AK
       CPG=AK*RG/(AK-1)
       TAUF=RHOP*DP**2/18.0
       SHR=CPS/CPG
       RLC=CPS/HLAT
       RADF=6*PLANK*EPSILON/(RHOP*CPS*DP)
       IPHASE=0
       IF(TP0.LT.TMELT) PS=0.0
       IF(TP0.GT.TMELT) PS=1.0
       DEN0=P0/RG/T0
       IC(1)=0
C---------------------------------------------------------------------
C  CALL SUBROUTINE FOR DUCT GEOMETRY
C
       CALL DUCTGEO(AL,DX,DIA1,DIA2,IMAX,XL,DIA,AR)
C---------------------------------------------------------------------
C--INITIALIZE CONSERVATIVE VARIABLES
       X0=DEN0*AR(1)*U0
       Y0=X0*U0+P0*AR(1)
       Z0=X0*(CPG*T0+U0*U0/2)
       DEN1=DEN0
       U1=U0
       T1=T0
       P1=P0
       X1=X0
       Y1=Y0
       Z1=Z0
       UP1=UP0
       TP1=TP0
       U(1)=U0
       T(1)=T0
       P(1)=P0
       UP(1)=UP0
       TP(1)=TP0
       IC(1)=0
C---------------------------------------------------------------------
C      BEGIN DO LOOP
C---------------------------------------------------------------------
       DO 100 ISTEP=2,IMAX
C--ESTIMATE VELOCITY, PRESSURE, DENSITY, AND TEMPERATURE AT NEXT STATION
       U2=U1*AR(ISTEP)/AR(ISTEP-1)
       P2=P1-DEN1*(U2*U2-U1*U1)/2
       DEN2=DEN1
       T2=T1
       IC(ISTEP)=0
       XL(ISTEP)=XL(ISTEP-1)+DX
C---------------------------------------------------------------------
C--BEGINNING ITERATION FOR PRESSURE AND VELOCITY AT NEXT STATION
C---------------------------------------------------------------------
   20  CONTINUE
       IC(ISTEP)=IC(ISTEP)+1
       U2E=U2
       UAVE=(U2+U1)*0.5
       DENAV=(DEN1+DEN2)*0.5
       TAVE=(T1+T2)*0.5
       VISC=VISCO*SQRT(TAVE/T0)
       DIAAV=(DIA(ISTEP)+DIA(ISTEP-1))*0.5
       RE=DENAV*UAVE*DIAAV/VISC
       RR=ROUGH/DIAAV
       FF=0.25/(0.434*ALOG(RR/3.7)+5.74/RE**0.9)**2
C---------------------------------------------------------------------
C--CALL THE PARTICLE SUBROUTINE TO EVALUATE PARTICLE VELOCITY AND
C  TEMPERATURE AND THE SOURCE TERMS
       CALL PARTICLE(UP1,TP1,U2,T2,UP2,TP2,PS,SMOMP,SENERP,IPHASE)
C---------------------------------------------------------------------
C      EVALUATING SOURCE TERMS
C---------------------------------------------------------------------
C
C--MASS SOURCE
       SMASS=0.0
C--MOMENTUM SOURCE
       SMOM=0.0
C------SOURCE TERM FOR AREA CHANGE
       SMOM=SMOM+(P1+P2)*(AR(ISTEP)-AR(ISTEP-1))*0.5
C------SOURCE TERM FOR FRICTION
       SMOM=SMOM-FF*DX*DENAV*UAVE*UAVE*3.1416*DIAAV
C------MOMENTUM SOURCE TERM FOR PARTICLES
       SMOM=SMOM+SMOMP
C--ENERGY SOURCE
       SENER=0.0
C------ENERGY SOURCE TERM FOR PARTICLES
       SENER=SENER+SENERP
C---------------------------------------------------------------------
C      UPDATING SOURCE TERMS
C---------------------------------------------------------------------
       X2=X1+SMASS
       Y2=Y1+SMOM
       Z2=Z1+SENER
       SONIC=FAC2*X2*Z2/Y2/Y2
       WRITE(*,2004) X2,Y2,Z2,SONIC,ISTEP
       IF(SONIC.GT.1.0) THEN
          WRITE(*,1001)
          STOP
       END IF
C---------------------------------------------------------------------
C      EVALUATING PRIMITIVE VARIABLES AT NEXT STATION
C---------------------------------------------------------------------
       U2=FAC1*Y2*(1.0-SQRT(1.0-SONIC))/X2
       DEN2=X2/AR(ISTEP)/U2
       P2=(Y2-X2*U2)/AR(ISTEP)
       T2=(Z2/X2-U2*U2/2)/CPG
       IF(IC(ISTEP).GT.20) THEN
            WRITE(*,1002)
            STOP
       END IF
       IF(ABS(U2E-U2).GT.RES) GOTO 20
C---------------------------------------------------------------------
C      END OF ITERATION LOOP AND INITIALIZING FOR NEXT STEP
C---------------------------------------------------------------------
       U1=U2
C       U(ISTEP)=U1
       U(ISTEP)=U2
       T1=T2
       T(ISTEP)=T2
       UP1=UP2
       UP(ISTEP)=UP2
       TP1=TP2
       TP(ISTEP)=TP2
       P1=P2
       P(ISTEP)=P2
       DEN1=DEN2
       DEN(ISTEP)=DEN2
C--CHECKING ON CHANGE OF PHASE
       IF(PS.EQ.0.0.AND.TP2.GT.TMELT) IPHASE=1
       IF(PS.EQ.1.0.AND.TP2.LT.TMELT) IPHASE=1
       IF(PS.GT.1.0) THEN
                PS=1.0
                IPHASE=0
       ENDIF
       IF(PS.LT.0.0) THEN
                PS=0.0
                IPHASE=0
       ENDIF
C---------------------------------------------------------------------
       X1=X2
       Y1=Y2
       Z1=Z2
  100  CONTINUE
C---------------------------------------------------------------------
C      PRINTOUT OF FINAL VALUES
C---------------------------------------------------------------------
       WRITE(1,2000) P0,T0,U0
       WRITE(1,2002)
       DO 200 I=1,IMAX,ISKIP
       WRITE(1,2003) XL(I),U(I),T(I),P(I),UP(I),TP(I),IC(I)
  200  CONTINUE
       STOP
 1001  FORMAT(' INITIAL VELOCITY TOO HIGH')
 1002  FORMAT(' TOO MANY ITERATIONS')
 2000  FORMAT(' GAS FLOW IN A DUCT USING CONSERVATIVE VARIABLES'
     1 /' PRES=',E12.4,4X,'TEMP=',F10.2,4X,'VEL=',F8.4)
 2002  FORMAT(/'    DISTANCE    VELOCITY   GAS TEMP   PRESSURE   PART
     1 VEL  PART TEMP  ITER')
 2003  FORMAT(3F12.4,1E12.4,2F12.4,I3)
 2004  FORMAT(4E12.4,I5)

       END
C---------------------------------------------------------------------
C
C      SUBROUTINE FOR DUCT GEOMETRY
C
C      THIS SUBROUTINE CREATES A COSINE CURVE FOR THE NOZZLE GEOMETRY
C
C---------------------------------------------------------------------
      SUBROUTINE DUCTGEO(AL,DX,DIA1,DIA2,IT,XL,DIAX,ARX)
      DIMENSION XL(IT),DIAX(IT),ARX(IT)
      DIAM=(DIA1+DIA2)*0.5
      DIAA=(DIA1-DIA2)*0.5
      PIL=2*3.1416/AL
      DO 10 I=1,IT
      XL(I)=(I-1)*DX
      DIAX(I)=DIAM+DIAA*COS(XL(I)*PIL)
      ARX(I)=0.785*DIAX(I)**2
   10 CONTINUE
      RETURN
      END
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C
C     SUBROUTINE FOR CALCULATING PARTICLE VELOCITY, TEMPERATURE AND
C     COUPLING TERMS
C
      SUBROUTINE PARTICLE(UP1,TP1,U2,T2,UP2,TP2,PS,SMOMP,SENERP,IPHASE)
C
C--CALCULATE VELOCITY USING THE IMPLICIT QUADRATIC APPROACH. THIS SCHEME
C     WILL NOT WORK IF VELOCITY GOES NEGATIVE!
C
      COMMON TAUF,VISC,AG,FRICP,DX,DIAAV,PRN,SHR,UAVE,DENAV,DP,TMELT,
     1       TW,RLC,RADF,X1,CPS,TAVE,ZL
C---------------------------------------------------------------------
C
      IF(UP1.LT.0.0) THEN
         WRITE(*,1003)
         STOP
      ENDIF
      RER=ABS(UAVE-UP1)*DENAV*DP/VISC
      IF(RER.LT.0.1) THEN
         DFAC=1.0
      ELSE
         DFAC=1.0+0.15*RER**0.67
      END IF
      TAU=TAUF/VISC
      FAC1=DFAC*DX/TAU
      FAC2=1+FRICP*DX*2/DIAAV
      UP2=(-FAC1+SQRT(FAC1**2+2*FAC2*(UP1**2/2+AG*DX+FAC1*U2)))/FAC2
C
C---------------------------------------------------------------------
C
      UPAVE=(UP1+UP2)*0.5
      IF(REF.LT.0.1) THEN
        ANU=2.0
      ELSE
        ANU=2.0+0.6*RER**0.5*PRN**0.33
      END IF
      FAC3=ANU*DX/(3*PRN*SHR*TAU*UPAVE)
      FAC4=RADF*DX/UPAVE
      IF(IPHASE.EQ.1) THEN
        TP2=TMELT
        PS=PS+FAC3*RLC*(T2-TMELT)-FAC4*RLC*(TP2**4-TW**4)
      ELSE
        TP2=(TP1+FAC3*T2)/(1.0+FAC3)
      END IF
      TPAVE=(TP1+TP2)*0.5
C-----EVALUATION OF THE SOURCE TERM FOR MOMENTUM
      SMOMP=ZL*X1*FAC1*(UPAVE-UAVE)/UPAVE
C-----EVALUATION OF SOURCE TERM FOR ENERGY - HEAT TRANSFER AND WORK
      SENERP=ZL*X1*FAC3*(TPAVE-TAVE)*CPS+UPAVE*SMOMP
C
      RETURN
 1003 FORMAT(' VELOCITY IS NEGATIVE-PARTICLE SUBROUTINE INOPERATIVE')
      END
C---------------------------------------------------------------------
