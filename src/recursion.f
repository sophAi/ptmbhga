        SUBROUTINE RECURSION(XINIB_ORG1,XINIB_MIN,NDIM,ESYS,EMINGL,
     &SEED,RATE,IBIN,BETA,AL,KEN)
        IMPLICIT NONE
        INTEGER  IBIN,NDIM,SEED,MIN_S,MIN_T,JMAX,JMAX_ORG,KBIN
        INTEGER KEN(IBIN),N0,MCMAX,STEPMAX,J1,J2,ERANGE_SW,SW_SAME_TEMP
        INTEGER DG_SW,BH_SW,view_sweep
        INTEGER SW_ENSNUM,SW_PARNUM,SW_SAME,SW_GA,NP(100),NO1,NO2
        INTEGER sw_selec_par,KCOUNT1,KCOUNT2,COUNTER1,COUNTER2,PRE_BH
        REAL*8 XINIB_ORG1(NDIM),XINIB_ORG2(NDIM),RND,RATE,DG,lmin,RMS,
     &FIT(100),FITG(101),RND_DIFF
        REAL*8 VMIN,VMAX,DMAX,VMIN_ORG,VMAX_ORG,DMAX_ORG,XINIB_MIN(NDIM)
        REAL*8 GRAD(NDIM),BETA(IBIN),AL(IBIN),EMINGL,
     &XINIB_ENS(100,NDIM)
        REAL*8 XINIB_NEW1(NDIM),XINIB_NEW2(NDIM),VT(2000),VTR(2000),
     &VT_ORG(2000),VTR_ORG(2000)
        REAL*8 EMAX,EMIN,DE,BETA0,BETANW,ALNW,BETAOD,ALOD,ESYS,ESYSNW1
        REAL*8 DIFF,SCALEFAC,TEMP,EMAX_TEMP,EMIN_TEMP,TDE,ESYSNW2,
     &POT_ENS(100),ESYSNW_TEMP,XINIB_TEMP(NDIM)
        LOGICAL EVAP,CFLAG
        COMMON/MULTIC1/EMAX,EMIN
        COMMON/MULTIC2/ERANGE_SW,PRE_BH
        COMMON/MULTIC3/DE
C        COMMON/GLOBAL/EMINGL
        COMMON/UPPER/ BETA0, N0,STEPMAX,MCMAX
        COMMON/TEMP1/ SCALEFAC,TEMP
        COMMON/vt1/ VT,VTR
        COMMON/vt2/ JMAX
        COMMON/vt3/ VMIN,VMAX,DMAX
        COMMON/guide1/ DG,DG_SW
        COMMON/BH1/ BH_SW
        COMMON/BH4/ TDE
        COMMON/SW_GA0/ SW_ENSNUM,SW_PARNUM,SW_SAME
        COMMON/SW_GA2/ SW_GA
        COMMON/VIEW1/ view_sweep
        RATE=0.D0
        MIN_S=0
        MIN_T=0
        J2=0
        KCOUNT2=0
        COUNTER1=SW_ENSNUM+1
C        view_sweep=1
        IF(BH_SW.NE.2.OR.BH_SW.NE.3)THEN
          EMAX_TEMP=EMAX
          EMIN_TEMP=EMIN
          DE=(EMAX_TEMP-EMIN_TEMP)/DBLE(IBIN)
          IF(ERANGE_SW.EQ.2) THEN 
            EMIN=EMINGL
          ELSE IF(ERANGE_SW.EQ.1) THEN 
            EMIN=EMINGL-(TDE/2.D0)
            EMAX=EMIN+TDE
          ENDIF
        ENDIF
        DO J1=1,NDIM
          XINIB_ORG2(J1)=XINIB_MIN(J1)
C          write(*,*) J1,"original xinib=",xinib_org2(J1),NDIM
        ENDDO
        EVAP=.FALSE.
        IF (STEPMAX.EQ.0) RETURN
        DO 9999 KCOUNT1=1,STEPMAX                ! Chose atoms NT times
C        N=RAN(SEED)*FLOAT(NDIM)+1              ! This chose an atom
C =======================NOW MOVE AN ATOM========================
C         XINIB(KCOUNT)=XINIB(KCOUNT)+(XINIB(KCOUNT)*RND(SEED)-1)
C =======================NOW CHECK IF THIS IS OK WITH B CONDITIONS ===
C        R2=XNEW*XNEW+YNEW*YNEW+ZNEW*ZNEW
C        IF(R2.GT.R2MAX) GOTO 998                  ! R2MAX MUST BE SET IN INIPAR
C =======================PASSED THE BC TEST, NOW EVALUATE ENERGY OF SYSTEM
          JMAX_ORG=JMAX
          VMIN_ORG=VMIN
          VMAX_ORG=VMAX
          DMAX_ORG=DMAX
          DO J1=1,NDIM
            VT_ORG(J1)=VT(J1)
            VTR_ORG(J1)=VTR(J1)
          ENDDO
C          DO J1=1,NDIM
C            XINIB_NEW1(J1)=XINIB_ORG1(J1)
C          ENDDO
C          GOTO 997
996       CALL step_move(XINIB_ORG1,XINIB_ORG2,XINIB_NEW1,XINIB_NEW2,
     &NDIM,SEED)   !use move to selet move method,ex mutation
997       CALL LOCAL_MIN(XINIB_NEW1,ESYSNW1,NDIM,SEED,CFLAG,MIN_S,
     &MIN_T,EVAP)
          IF (EVAP) THEN
            JMAX=JMAX_ORG
            VMIN=VMIN_ORG
            VMAX=VMAX_ORG
            DMAX=DMAX_ORG
            DO J1=1,NDIM
              VT(J1)=VT_ORG(J1)
              VTR(J1)=VTR_ORG(J1)
            ENDDO
            J2=J2+1
            IF(J2.ge.5) GOTO 998
            GOTO 996
          ENDIF
          ESYSNW2=0.D0
          IF(SW_GA.eq.1)then
            CALL LOCAL_MIN(XINIB_NEW2,ESYSNW2,NDIM,SEED,CFLAG,MIN_S,
     &MIN_T,EVAP)
C          write(*,*) ESYSNW,ESYS,GBIN_TOTAL
            IF (EVAP) THEN
              JMAX=JMAX_ORG
              VMIN=VMIN_ORG
              VMAX=VMAX_ORG
              DMAX=DMAX_ORG
              DO J1=1,NDIM
                VT(J1)=VT_ORG(J1)
                VTR(J1)=VTR_ORG(J1)
              ENDDO
              J2=J2+1
              IF(J2.ge.5) GOTO 998
              GOTO 996
            ENDIF
          ENDIF
          IF (ESYSNW2.LT.ESYSNW1)THEN
            ESYSNW_TEMP=ESYSNW2
            ESYSNW2=ESYSNW1
            ESYSNW1=ESYSNW_TEMP
            DO J1=1,NDIM
              XINIB_TEMP(J1)=XINIB_NEW2(J1)
              XINIB_NEW2(J1)=XINIB_NEW1(J1)
              XINIB_NEW1(J1)=XINIB_TEMP(J1)
            ENDDO
          ENDIF
C=============START THE METROPOLIS PART=================
C=============START MMC=================================
          IF(BH_SW.EQ.1)THEN
            IF (ESYSNW1.GT.EMAX_TEMP)THEN
              BETANW=BETA0
              ALNW  =0.D0
            ELSE IF (ESYSNW1.LE.EMIN_TEMP)THEN
              BETANW=BETA(IBIN)
              ALNW=AL  (IBIN)
            ELSE
              KBIN=INT((EMAX_TEMP-ESYSNW1)/DE+1.D0)
              BETANW=BETA(KBIN)
              ALNW=AL  (KBIN)
            ENDIF
            IF (ESYS.GT.EMAX_TEMP)THEN
              BETAOD=BETA0
              ALOD  =0.D0
            ELSE IF (ESYS.LE.EMIN_TEMP)THEN
              BETAOD=BETA(IBIN)
              ALOD=AL  (IBIN)
            ELSE
              KBIN=INT((EMAX_TEMP-ESYS)/DE+1.D0)
              BETAOD=BETA(KBIN)
              ALOD=AL  (KBIN)
            ENDIF
            DIFF=-BETANW*ESYSNW1-ALNW+BETAOD*ESYS+ALOD
C=============END OF MMC================================
C=============START OF AMC==============================
          ELSE IF(BH_SW.EQ.5)THEN
            IF (ESYSNW1.GT.EMAX_TEMP)THEN
              DIFF=-999999.D0
            ELSE IF (ESYSNW1.LE.EMIN_TEMP)THEN
              KBIN=INT(EMAX_TEMP-EMINGL/DE+1.D0)
              DIFF=-AL(KBIN)-BETA(IBIN)*(ESYSNW1-EMINGL)
            ELSE
              KBIN=INT(EMAX_TEMP-ESYSNW1/DE+1.D0)
              DIFF=-AL(KBIN)
            ENDIF
C==============END OF AMC===============================
C==============START OF PGA=============================
          ELSE IF(BH_SW.EQ.3)THEN
            DIFF=1.D0
            IF (KCOUNT2.LE.SW_ENSNUM)THEN
C              WRITE(*,*) KCOUNT2
              KCOUNT2=KCOUNT2+1
              DO J1=1,NDIM
                XINIB_ENS(KCOUNT2,J1)=XINIB_NEW1(J1)
C                write(*,*) KCOUNT2,J1,"XINIB=",XINIB_ENS(KCOUNT2,J1)
                IF(SW_GA.EQ.1.AND.(KCOUNT2+1).LE.SW_ENSNUM)THEN
                  XINIB_ENS(KCOUNT2+1,J1)=XINIB_NEW2(J1)
C               write(*,*) KCOUNT2+1,J1,"XINIB=",XINIB_ENS(KCOUNT2+1,J1)
                ENDIF
              ENDDO
              POT_ENS(KCOUNT2)=ESYSNW1
              IF(SW_GA.EQ.1)THEN
                KCOUNT2=KCOUNT2+1
                POT_ENS(KCOUNT2)=ESYSNW2
              ENDIF
C             Collect 20 ensemble!
            ELSE
              IF (COUNTER1.EQ.SW_ENSNUM+1)THEN
                CALL sw_ene_sort(POT_ENS,SW_SAME_TEMP,SW_ENSNUM,NP)
              ENDIF
              COUNTER1=COUNTER1-1
              DO J1=1,NDIM
                XINIB_ENS(NP(COUNTER1),J1)=XINIB_NEW1(J1)
              ENDDO
              POT_ENS(NP(COUNTER1))=ESYSNW1
              IF(SW_GA.EQ.1)THEN
                COUNTER1=COUNTER1-1
                DO J1=1,NDIM
                  XINIB_ENS(NP(COUNTER1),J1)=XINIB_NEW2(J1)
                ENDDO
                POT_ENS(NP(COUNTER1))=ESYSNW2
              ENDIF
C              WRITE(*,*) "SWGA=",SW_GA
C              WRITE(*,*) COUNTER1,ESYSNW1,ESYSNW2,POT_ENS(NP(COUNTER1))
              DO J1=1,SW_ENSNUM
C                WRITE(*,*) J1,",POT=",POT_ENS(NP(J1)),COUNTER1,KCOUNT2
              ENDDO
C              PAUSE
              IF(COUNTER1.LE.SW_PARNUM+1)THEN
                CALL sw_ene_sort(POT_ENS,SW_SAME_TEMP,SW_ENSNUM,NP)
                COUNTER1=SW_ENSNUM+1
                IF (SW_SAME_TEMP.GE.SW_SAME) THEN
                  WRITE(*,*) "SAME=",SW_SAME_TEMP,",EMINGL=",EMINGL
                  RETURN
                ENDIF
              ENDIF 
              CALL sw_fitness(POT_ENS,SW_PARNUM,SW_ENSNUM,FIT,FITG,NP)
              NO1=NP(sw_selec_par(SW_PARNUM,SW_ENSNUM,FITG,NP,SEED))
              NO2=NP(sw_selec_par(SW_PARNUM,SW_ENSNUM,FITG,NP,SEED))
              DO J1=1,NDIM
                XINIB_ORG1(J1)=XINIB_ENS(NO1,J1)
                XINIB_ORG2(J1)=XINIB_ENS(NO2,J1)
              ENDDO
            ENDIF               
          ELSE
            DIFF=((-ESYSNW1+ESYS)/TEMP)
          ENDIF
          IF (DIFF.GT.0.D0)GOTO 999                        ! ACCEPTING
          RND_DIFF=RND(SEED)
          IF (DEXP(DIFF).GT.RND_DIFF)GOTO 999             ! ACCEPTING
C          WRITE(*,*) "REJECT,",KCOUNT,ESYS,ESYSNW,TEMP
          JMAX=JMAX_ORG
          VMIN=VMIN_ORG
          VMAX=VMAX_ORG
          DMAX=DMAX_ORG
          DO J1=1,NDIM
            VT(J1)=VT_ORG(J1)
            VTR(J1)=VTR_ORG(J1)
          ENDDO
          if(view_sweep.eq.1)then
      write(*,*) "BH STEP=",KCOUNT1,",Local min=",ESYSNW1,
     &" ,REJECT! DE=",DEXP(DIFF),"<",RND_DIFF
          endif
          GOTO 998                                        ! REJECT
C========================================================
999       CONTINUE
C==================ACCEPTING THE MOVE====================
C          WRITE(*,*) "ACCEPT IN SEEEP"
C          WRITE(*,*) "ACCEPT,",KCOUNT,ESYS,ESYSNW,TEMP
          ESYS=ESYSNW1
          DO J1=1,NDIM
            XINIB_ORG1(J1)=XINIB_NEW1(J1)
          ENDDO
          RATE=RATE+1.D0
          if(view_sweep.eq.1)then
            write(*,*) "BH STEP=",KCOUNT1,",Local min=",ESYSNW1,
     &" ,ACCEPT! HIT RATE=",RATE/DBLE(STEPMAX)
          endif
C999
C998       CONTINUE
C========ACCEPT OR REJECT...NOW STATISTICS===============
C FINDING NEW MINIMUM
998       IF (ESYS.LT.EMINGL) THEN
C            WRITE(*,*) "ESYS=",ESYS,",EMINGL=",EMINGL
            EMINGL=ESYS
C            write(*,*) ESYS,KCOUNT
            DO J1=1,NDIM
              XINIB_MIN(J1)=XINIB_ORG1(J1)
              XINIB_ORG2(J1)=XINIB_MIN(J1)
            ENDDO
          ENDIF
C ADDING UP THE HISTORGR/
          IF (ESYS.LT.EMAX_TEMP.AND.ESYS.GT.EMIN_TEMP)THEN
            KBIN=int((EMAX_TEMP-ESYS)/DE+1.D0)
c        print *, kbin, emax,esys
c        pause
            KEN(KBIN)=KEN(KBIN)+1
C            write(*,*) KBIN,".KEN=",KEN(KBIN)
          ENDIF
          IF (BH_SW.NE.2.OR.BH_SW.NE.3)THEN
            IF (ERANGE_SW.eq.2)THEN
              IF(EMINGL.LT.EMIN)THEN
                EMIN=EMINGL
              ENDIF
              EMAX=dble(max(ESYSNW1,(EMIN+TDE)))
            ELSE IF (ERANGE_SW.eq.1)THEN
              IF(EMINGL.LT.EMIN)THEN
                EMIN=(EMINGL-(TDE/2.D0))
              ENDIF
              EMAX=dble(max(ESYSNW1,(EMIN+TDE)))
            ENDIF
          ENDIF          
9999    CONTINUE
C       PRINT *, ESYS,
        RATE=RATE/DBLE(STEPMAX)
C        call POTENTIAL(ndim,xinib_min,grad,lmin,RMS,.false.,EVAP)
C        if (dabs(lmin-EMINGL).gt.0.00001D0)then
C        write(*,*) "final output local min in sweep=",lmin,",EMINGL=",
C     &EMINGL
C         pause "ERRPR!!!!"
C        endif
C        DO J1=1,NDIM
C          WRITE(*,*) J1,XINIB_MIN(J1)
C        ENDDO
C        pause
        RETURN
        END

C        SUBROUTINE VMDFMT
C        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        PARAMETER (IBIN=20 )                      ! No. of multicanonical bins
C        COMMON/M_PAR/R2MAX, ZETA0, CETA0, P, Q
C        COMMON/MULTIC/EMAX,EMIN,DE,ESYS
C        COMMON/UPPER /BETA0, N0,MCMAX,MULMAX
C        character *6, conect
C        conect='CONECT'
C        do 5 J=1,NT
C        write(51,'(a14)')'HEADER    PROTEIN'
C        write(51,33)'ATOM', J,     'CA', 'UNK', J,x(j),y(j),z(j)
C 5      continue
C 33     FORMAT      (a4,I7, 2X,a2,2X,  a3, I6,4X,3(F8.3))
C        write(51,34)'END   '
C 34     FORMAT(a6,3(I5))
C        end

*****

        SUBROUTINE NEWAB(IBIN,EMINGL,KEN,BETA,AL,S,G_ALL)
        IMPLICIT NONE
        INTEGER IBIN,N0,STEPMAX,MCMAX,I,IMIN,KBIN
        INTEGER KEN(IBIN),BH_SW,G_SW
        REAL*8 BETA(IBIN),AL(IBIN),S(IBIN),BETA0,EMINGL,DE,EMAX,EMIN
     &,EPRE,EP,BETDIF,G_ALL(IBIN),G(IBIN)
C        CHARACTER name*6
        COMMON/MULTIC1/EMAX,EMIN
        COMMON/MULTIC3/DE
C        COMMON/GLOBAL/EMINGL
        COMMON/UPPER/ BETA0, N0,STEPMAX,MCMAX
        COMMON/BH1/ BH_SW
        COMMON/WEIGHT1/G_SW
C        COMMON/name1/ name
c       BETA0=BETA0+DBETA
        DE=(EMAX-EMIN)/DBLE(IBIN)
        DO I=1,IBIN
          IF (KEN(I).GE.20) THEN
            S(I)=S(I)+ DLOG(DBLE(KEN(I)))
          ENDIF
          IF (G_SW.EQ.1)THEN
            IF (I.GT.1.AND.KEN(I-1).NE.0.AND.KEN(I).NE.0)THEN
              G(I)=DBLE(KEN(I-1)*KEN(I))/DBLE(KEN(I-1)+KEN(I))
              G_ALL(I)=G(I)+G_ALL(I)
              G(I)=G(I)/G_ALL(I)
            ELSE
              G(I)=0.D0
            ENDIF
          ELSE
            G(I)=1.D0
          ENDIF
        ENDDO
        IF(BH_SW.EQ.1)THEN
          BETA(1)=BETA0
          AL(1)  =0.D0 
          IF(EMINGL.LE.EMAX.AND.EMINGL.GE.EMIN)THEN 
            IMIN=int((EMAX-EMINGL)/DE+1.D0)
            DO  I=2,IMIN
              BETA(I)=BETA0 + ((S(I-1)-S(I))/DE)*G(I-1)
            ENDDO
            BETDIF=BETA(IMIN)-BETA(IMIN-1)
            DO I=IMIN+1,IBIN
              BETA(I)=BETA(I-1)
ccc+BETDIF
            ENDDO
            DO I=2,IBIN
              EPRE=EMAX-DBLE(I-2)*DE
              AL(I)=AL(I-1)+ (BETA(I-1)-BETA(I))*EPRE
            ENDDO
          ELSE
            WRITE(*,*) "WARNING!!!!!!!...LOCAL MINIMUM IS OUT OF RANGE"
            WRITE(*,*) "YOU HAVE TO RESELECT THE ENERGY RANGE"
            WRITE(*,*) "THE MULTI-CANONICAL BH WILL FAIL!!!"
            WRITE(*,*) "EMIN=",EMIN,",lmin=",EMINGL,",EMAX=",EMAX
            WRITE(*,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
          ENDIF
        ENDIF
        IF(BH_SW.EQ.5)THEN
          DO I=1,IBIN
            BETA(I)=0.D0
          ENDDO
          IMIN=int((EMAX-EMINGL)/DE+1.D0)
          BETA(IBIN)=(S(1)-S(IMIN))/(EMAX-EMINGL)
        ENDIF
C        PRINT *, 'Emin=',EMINGL, 'bin no=',IMIN

C        PRINT *, ' I    BETA(I)   ALPHA(I)  KEN(I)    E(I) '
C        DO  I=1,IBIN
C          EP=EMAX-DBLE(I-1)*DE
C        ENDDO
C60        WRITE(6,733) i, beta(i), al(i), ken(i),EP
C 733    FORMAT(I6,F9.4,1X, F9.4, 1X, I8, F9.4)
C        pause
C        OPEN(46,file=name//".thm",status="replace")
C        DO 71 I=1,NDIM
C 71     WRITE(46,*)XINIB(I)

C        do 70 i=1,ibin
C 70     write(46,*) i, beta(i), al(i), s(i)
C        write(46,*) esys, ' =Esystem'
C        close(46)
C        pause
        DO  I=1,IBIN
          KEN(I)=0
        ENDDO
        RETURN
        END

