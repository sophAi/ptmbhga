        SUBROUTINE PROTEIN(NATOMS,X,GRAD,ENR,VT,EVAP,GRADT)
        INTEGER NATOMS,J1,j2,J3,j4,SW,I_TYPE(10000),TYPE
        REAL*8 VT(NATOMS),DIST,DIST1,DIST2,DIST3,CUTOFF,GRAD(NATOMS*3)
        REAL*8 FACTOR1,FACTOR2,RADIUS,NONE,ECOIL,G(NATOMS,NATOMS)
        REAL*8 COS_THEB,THEB,THETA0,CUTOFF_LENGTH,MAG_BA,MAG_BC,MAG_ABC
        REAL*8 X(NATOMS*3),ENR_LJ,ENR_BOND1,R3,A0,A1,A2,A3,A4,A5,B1,B2
        REAL*8 D11,D12,D13,D21,D22,D23,D31,D32,D33,SQR_ABC,ENR_DIH,ENR
        REAL*8 CROSS_ABCX,CROSS_ABCY,CROSS_ABCZ,MAG_CRO_AB
        REAL*8 CROSS_BCDX,CROSS_BCDY,CROSS_BCDZ,MAG_CRO_BC
        REAL*8 CROSS_ABCDX,CROSS_ABCDY,CROSS_ABCDZ,MAG_ABCD
        REAL*8 CRO_AB_DOT_CRO_BC,SIN_PHID,SIN_PHID2,SIN_PHID3
        REAL*8 COS_PHID,COS_PHID2,COS_PHID3,COS3DTH
        REAL*8 ACOIL,BCOIL,CCOIL,DCOIL,DSQ2,ENR_BOND2
        REAL*8 PARTIAL_IX,PARTIAL_IY,PARTIAL_IZ
        REAL*8 PARTIAL_JX,PARTIAL_JY,PARTIAL_JZ
        REAL*8 PARTIAL_KX,PARTIAL_KY,PARTIAL_KZ
        LOGICAL EVAP,GRADT
        COMMON /cgsw/ SW
        COMMON /engtype/ TYPE
        COMMON /INTER_TYPE1/ I_TYPE
        COMMON /LONG_RANGE1/ CUTOFF,ECOIL
        COMMON /accuracy3/ RADIUS,NONE
        COMMON /BOND1/ THETA0,CUTOFF_LENGTH
        COMMON /DIH1/ ACOIL,BCOIL,CCOIL,DCOIL,R3,A0,A1,A2,A3,A4,A5,B1,B2
C   I_TYPE=1   B
C   I_TYPE=2   L
C   I_TYPE=3   N
****************<<<LJ-E>>>************************
        ENR_LJ=0.D0
        EVAP=.FALSE.
        DO J1=1,NATOMS
          J3=J1*3
          VT(J1)=0.0D0
          DO J2=J1+1,NATOMS
            J4=J2*3
            GRAD(J3-2)=0.D0
            GRAD(J3-1)=0.D0
            GRAD(J3)=0.D0
            DIST1=(X(J3)-X(J4))**2.
            DIST2=(X(J3-1)-X(J4-1))**2.
            DIST3=(X(J3-2)-X(J4-2))**2.
       IF(DIST1.LT.CUTOFF.AND.DIST2.LT.CUTOFF.AND.DIST3.LT.CUTOFF)THEN
              DIST=DIST1+DIST2+DIST3
              IF (DSQRT(DIST).GT.RADIUS)THEN
                EVAP=.TRUE.
                IF (SW.EQ.11)RETURN
                ENR_LJ=ENR_LJ+(DIST-RADIUS)**2.
              ENDIF
              IF (I_TYPE(J1).EQ.1.AND.I_TYPE(J2).EQ.1)THEN
C  BB
                FACTOR1=1.D0
                FACTOR2=1.D0
                ELSE IF (I_TYPE(J1).EQ.2.OR.I_TYPE(J2).EQ.2)THEN
C  LB,BL,LL       
                  IF (TYPE.EQ.3)THEN
                    FACTOR1=2.D0/3.D0
                    FACTOR2=(-2.D0)/3.D0
                  ENDIF
                  IF (TYPE.EQ.4)THEN
                    FACTOR1=1.D0
                    FACTOR2=0.D0
                  ENDIF
                ELSE
C  NB,BN,NL,LN,NN
                  IF (TYPE.EQ.3)THEN
                    FACTOR1=1.D0
                    FACTOR2=0.D0
                  ENDIF
                  IF (TYPE.EQ.4)THEN
                    FACTOR1=1.D0
                    FACTOR2=0.D0
                  ENDIF
                ENDIF
                ENR_LJ=ENR_LJ+4.D0*ECOIL*((FACTOR1*(DIST**(-6.)))-
     &(FACTOR2*(DIST**(-3.))))
                VT(J1)=VT(J1)+4.D0*ECOIL*((FACTOR1*(DIST**(-6.)))-
     &(FACTOR2*(DIST**(-3.))))
                IF (GRADT)THEN 
                  G(J2,J1)=4.D0*ECOIL*((-12.D0)*FACTOR1*(DIST**(-14.))+
     &(-6.D0)*FACTOR2*(DIST**(-8.)))
                  G(J1,J2)=G(J2,J1)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          
              
****************<<<BOND-ANGLE>>>************************
          ENR_BOND1=0.D0
          ENR_BOND2=0.D0
          DO J1=1,NATOMS-2
            J3=J1*3
            COS_THEB=(X(J3)-X(J3+3))*(X(J3+6)-X(J3+3))
     &         +(X(J3-1)-X(J3+2))*(X(J3+5)-X(J3+2))
     &         +(X(J3-2)-X(J3+1))*(X(J3+4)-X(J3+1))
            MAG_BA=(X(J3)-X(J3+3))**2.+(X(J3-1)-X(J3+2))**2.+
     &           (X(J3-2)-X(J3+1))**2.
            MAG_BC=(X(J3+6)-X(J3+3))**2.+(X(J3+5)-X(J3+2))**2.
     &          +(X(J3+4)-X(J3+1))**2.
            MAG_ABC=DSQRT(MAG_BA*MAG_BC)
            THEB=DACOS(COS_THEB/MAG_ABC)
            ENR_BOND1=ENR_BOND1+10.D0*(THEB-THETA0)**2.
C*****BOND DISTANCE AND ORDER CONSTRAIN*************
            DIST1=(X(J3+1)-X(J3-2))**2+(X(J3+2)-X(J3-1))**2+
     &(X(J3+3)-X(J3))**2
            ENR_BOND1=ENR_BOND1+5.D0*(DIST1-CUTOFF_LENGTH)**2.
            IF (J1.EQ.NATOMS)THEN
              DIST2=(X(J3+4)-X(J3+1))**2.+(X(J3+5)-X(J3+2))**2.+
     &(X(J3+6)-X(J3+3))**2.
              ENR_BOND1=ENR_BOND1+5.D0*(DIST2-CUTOFF_ANGLE)**2.
            ENDIF
C*************END,START GRADIENT***************************
            IF (GRADT)THEN
C    I_TYPE
                PARTIAL_IX=(X(J3+4)-X(J3+1))/MAG_ABC+(X(J3-2)-X(J3+1))*
     &(COS_THEB/MAG_ABC)/MAG_BA
                PARTIAL_IY=(X(J3+5)-X(J3+2))/MAG_ABC+(X(J3-1)-X(J3+2))*
     &(COS_THEB/MAG_ABC)/MAG_BA
                PARTIAL_IZ=(X(J3+6)-X(J3+3))/MAG_ABC+(X(J3)-X(J3+3))*
     &(COS_THEB/MAG_ABC)/MAG_BA           
C    J_TYPE
                PARTIAL_JX=((2.D0*X(J3+1)-X(J3-2)-X(J3+4))/MAG_ABC)-
     &((X(J3-2)-X(J3+1))*(COS_THEB/MAG_ABC)/MAG_BA)-
     &((X(J3+4)-X(J3+1))*(COS_THEB/MAG_ABC)/MAG_BC)
                PARTIAL_JY=(2.D0*X(J3+2)-X(J3-1)-X(J3+5))/MAG_ABC-
     &((X(J3-1)-X(J3+2))*(COS_THEB/MAG_ABC)/MAG_BA)-
     &((X(J3+5)-X(J3+2))*(COS_THEB/MAG_ABC)/MAG_BC)
                PARTIAL_JZ=(2.D0*X(J3+3)-X(J3)-X(J3+6))/MAG_ABC-
     &((X(J3)-X(J3+3))*(COS_THEB/MAG_ABC)/MAG_BA)-
     &((X(J3+6)-X(J3+3))*(COS_THEB/MAG_ABC)/MAG_BC)
C    K_TYPE
                PARTIAL_KX=(X(J3-2)-X(J3+1))/MAG_ABC+(X(J3+4)-X(J3+1))*
     &(COS_THEB/MAG_ABC)/MAG_BC
                PARTIAL_KY=(X(J3-1)-X(J3+2))/MAG_ABC+(X(J3+5)-X(J3+2))*
     &(COS_THEB/MAG_ABC)/MAG_BC
                PARTIAL_KZ=(X(J3)-X(J3+3))/MAG_ABC+(X(J3+6)-X(J3+3))*
     &(COS_THEB/MAG_ABC)/MAG_BC
                
                ENR_BOND2=(-2.D0*(THEB-THETA0)/DSIN(THEB))
                
                GRAD(J3-2)=GRAD(J3-2)+ENR_BOND2*PARTIAL_IX
                GRAD(J3-1)=GRAD(J3-1)+ENR_BOND2*PARTIAL_IY
                GRAD(J3)=GRAD(J3)+ENR_BOND2*PARTIAL_IZ
                GRAD(J3+1)=GRAD(J3+1)+ENR_BOND2*PARTIAL_JX
                GRAD(J3+2)=GRAD(J3+2)+ENR_BOND2*PARTIAL_JY
                GRAD(J3+3)=GRAD(J3+3)+ENR_BOND2*PARTIAL_JZ
                GRAD(J3+4)=GRAD(J3+4)+ENR_BOND2*PARTIAL_KX
                GRAD(J3+5)=GRAD(J3+5)+ENR_BOND2*PARTIAL_KY
                GRAD(J3+6)=GRAD(J3+6)+ENR_BOND2*PARTIAL_KZ
            ENDIF
          ENDDO


***************<<DIH_ANGLE>>>****************************
              
        ENR_DIH=0.D0
        DO J1=1,NATOMS-3
          J3=J1*3
          D11=X(J3+1)-X(J3-2)
          D12=X(J3+2)-X(J3-1)
          D13=X(J3+3)-X(J3)
          D21=X(J3+4)-X(J3+1)
          D22=X(J3+5)-X(J3+2)
          D23=X(J3+6)-X(J3+3)
          D31=X(J3+7)-X(J3+4)
          D32=X(J3+8)-X(J3+5)
          D33=X(J3+9)-X(J3+6)

          CROSS_ABCX=D12*D23-D13*D22
          <F6>CROSS_ABCY=D13*D21-D11*D23
          CROSS_ABCZ=D11*D22-D12*D21

          CROSS_BCDX=D33*D22-D32*D23
          CROSS_BCDY=D31*D23-D33*D21
          CROSS_BCDZ=D32*D21-D31*D22

          CROSS_ABCDX=CROSS_ABCY*CROSS_BCDZ-CROSS_ABCZ*CROSS_BCDY
          CROSS_ABCDY=CROSS_ABCZ*CROSS_BCDX-CROSS_ABCX*CROSS_BCDZ
          CROSS_ABCDZ=CROSS_ABCX*CROSS_BCDY-CROSS_ABCY*CROSS_BCDX

          MAG_ABCD=CROSS_ABCDX**2.+CROSS_ABCDY**2.+CROSS_ABCDZ**2.

          CRO_AB_DOT_CRO_BC=CROSS_ABCX*CROSS_BCDX
     &                     +CROSS_ABCY*CROSS_BCDY
     &                     +CROSS_ABCZ*CROSS_BCDZ

          MAG_CRO_AB=CROSS_ABCX**2.+CROSS_ABCY**2.+CROSS_ABCZ**2.

          MAG_CRO_BC=CROSS_BCDX**2.+CROSS_BCDY**2.+CROSS_BCDZ**2.

          SQR_ABC=DSQRT(MAG_CRO_AB*MAG_CRO_BC)

          SIN_PHID=DSQRT(MAG_ABCD)/SQR_ABC
          SIN_PHID2=SIN_PHID**2.
          SIN_PHID3=SIN_PHID2*SIN_PHID

          COS_PHID=(CRO_AB_DOT_CRO_BC)/SQR_ABC
          COS_PHID2=COS_PHID**2.
          COS_PHID3=COS_PHID*COS_PHID2

          COS3DTH=4.D0*COS_PHID3-3.D0*COS_PHID
          IF (I_TYPE(J1).NE.3.AND.I_TYPE(J1+1).NE.3.AND.
     &    I_TYPE(J1+2).NE.3)THEN
C***************Beta*************
            IF (TYPE.EQ.3)THEN
              ENR_DIH=ENR_DIH+ACOIL*(1.D0+COS_PHID)+BCOIL*(1.D0+COS3DTH)
            ENDIF
C***************Alpha************
            IF (TYPE.EQ.4)THEN
              DSQ2=DSQRT(2.D0)
              ENR_DIH=ENR_DIH+ACOIL*(1.D0-COS_PHID)+BCOIL*(1.D0+COS3DTH)
     &  +CCOIL*(1.D0+(1.D0/DSQ2)*COS_PHID-(1.D0/DSQ2)*SIN_PHID)
            ENDIF
          ENDIF
          IF ((I_TYPE(J1).EQ.3.AND.I_TYPE(J1+1).EQ.3.AND.
     & I_TYPE(J1+2).NE.3.AND.I_TYPE(J1+3).NE.3).OR.(I_TYPE(J1).EQ.
     & 3.AND.I_TYPE(J1+1).EQ.3.AND.I_TYPE(J1+2).EQ.3))THEN     
            IF (TYPE.EQ.3.OR.TYPE.EQ.4)THEN
              ENR_DIH=ENR_DIH+DCOIL*(1.D0+COS3DTH)
            ENDIF
          ENDIF
          IF (TYPE.EQ.5.AND.((I_TYPE(J1).EQ.3.AND.I_TYPE(J1+1).EQ.3).OR.
     &(I_TYPE(J1+2).EQ.3.AND.I_TYPE(J1+3).EQ.3)))THEN
            ENR_DIH=ENR_DIH+ACOIL*(1.D0+COS_PHID)+BCOIL*(1.D0+COS3DTH)
          ELSE
            ENR_DIH=ENR_DIH+A0+A1*COS_PHID+A2*COS_PHID2+A3*COS_PHID3
     &         +A4*(COS_PHID2**2.)+A5*COS_PHID2*COS_PHID3
     &         +B1*SIN_PHID3+B2*SIN_PHID3*SIN_PHID2
          ENDIF

C*****************<<<GRADIENT OF DIH>>>************************
          IF (GRADT)THEN
            GRAD_MAB=
            GRAD_MBC=
            GRAD_RAB_DOT_RBC=
            GRAD_RBC_DOT_RAB=
            GRAD_RAB_CRO_RBC=
            GRAD_RBC_CRO_RAB=
            GRAD_MAG_IABC=(-1.D0)*(((MAG_CRO_AB**2)*(MAG_CRO_BC**2))*
     &(DSQRT(MAG_CRO_BC)*GRAD_MAB+DSQRT(MAG_CRO_AB)*GRAD_MBC)
            GRAD_COS=CRO_AB_DOT_CRO_BC*GRAD_MAG_IABC+(1.D0/SQR_ABC)*
     &(GRAD_RAB_DOT_RBC+GRAD_RBC_DOT_RAB))
            IF (TYPE.EQ.3)THEN
              GRAD=(ACOIL-3.D0*BCOIL)*GRAD_COS+12.D0*BCOIL*COS_PHID2*
     &GRAD_COS
            ENDIF
            IF (TYPE.EQ.4)THEN
              GRAD_SIN=
              GRAD=((1.D0/DSQ2)*CCOIL-ACOIL-3.D0*BCOIL)*GRAD_COS+12.D0*
     &BCOIL*COS_PHID2*GRAD_COS+(1.D0/DSQ2)*CCOIL*GRAD_SIN
            ENDIF
          ENDIF 
       ENDDO
       ENR=ENR_LJ+ENR_BOND1+ENR_DIH
       
       RETURN
       END
