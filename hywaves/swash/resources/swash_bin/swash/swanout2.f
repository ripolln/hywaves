!     SWAN/OUTPUT       file 2 of 2
!
!  Contents of this file:
!     SWBLOK
!     SBLKPT
!     SWBLKP
!     SWTABP
!     SWRMAT
!
!***********************************************************************
!                                                                      *
      SUBROUTINE SWBLOK ( RTYPE, OQI , IVTYP, FAC, PSNAME, MXK, MYK,
     &                    IRQ  , VOQR, VOQ  , VOQK )                      40.51 40.31
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE SwashCommdata3, ONLY: KMAX
      USE SwashTimecomm
      USE outp_data                                                       40.13
!
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWASH team                               |
!   --|-----------------------------------------------------------|--
!
!
!     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
!     Copyright (C) 2010-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     34.01: Jeroen Adema
!     40.03: Nico Booij
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. UPDATE
!
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Nov. 99: NVAR in write statement replaced by OREQ(18)
!     40.13, Oct. 01: longer output filenames now obtained from array
!                     OUTP_FILES (in module OUTP_DATA)
!     40.30, May  03: extension to write block output to Matlab files
!     40.31, Jul. 03: small correction w.r.t. length of OVSNAM in
!                     call SWRMAT
!     40.31, Dec. 03: removing POOL construction
!     40.41, Jun. 04: some improvements with respect to MATLAB
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!       Preparing output in the form of a block that is printed by
!         subroutine SBLKPT
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       RTYPE   ch*4   input    type of output request:
!                               'BLKP' for output on paper,
!                               'BLKD' and 'BLKL' for output to datafile
!       PSNAME  ch*8   input    name of outpu frame
!       MXK     int    input    number of grid points in x-direction
!       MYK     int    input    number of grid points in y-direction
!       VOQR
!       VOQ
!
!  5. SUBROUTINES CALLING
!
!       SWOUTP (SWAN/OUTP)
!
!  6. SUBROUTINES USED
!
!       SBLKPT, SCUNIT, SFLFUN (all SWAN/OUTP), TABHED,
!       MSGERR, COPYCH, FOR
!       SWRMAT, TXPBLA
!
      LOGICAL STPNOW                                                      34.01
!
!  7. ERROR MESSAGES
!
!       If the point set is not of the type frame, an error message
!       is printed and control returns to subroutine OUTPUT
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       If output is on paper then
!           Call TABHED to print heading
!       Else
!           Call FOR to open file
!       ----------------------------------------------------------------
!       For each required variable do
!           Determine type of variable and factor of multiplication
!           Call SBLKPT to write block output to printer or datafile
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      CHARACTER (LEN=8) :: PSNAME       ! name of output locations        40.13
      CHARACTER (LEN=4) :: RTYPE        ! output type                     40.13
      INTEGER   VOQR(*), IPD
      INTEGER   OQI(4), IVTYP(OQI(3))                                     40.31
      REAL      VOQ(MXK*MYK,*), VOQK(MXK*MYK,0:KMAX,*), FAC(OQI(3))       40.31
      INTEGER IF, IL                                                      40.41 40.30
      INTEGER, SAVE :: IREC(MAX_OUTP_REQ)=0                               40.51 40.41
      LOGICAL, SAVE :: MATLAB=.FALSE.                                     40.41 40.30
      LOGICAL IVSTAT1, IVSTAT2, IVSTAT3
      CHARACTER (LEN=20) :: CTIM                                          40.41
      CHARACTER (LEN=30) :: NAMVAR                                        40.41
      CHARACTER (LEN=4)  :: KNUM                                          SWASH
      CHARACTER (LEN=6)  :: CRFORM

      INTEGER, SAVE :: IENT=0                                             40.13
      IF (LTRACE) CALL STRACE (IENT,'SWBLOK')
!
!     **** obtain destination and number of variables from array OUTR ***
      NREF = OQI(1)                                                       40.31
      IF (RTYPE .EQ. 'BLKP') THEN
!       printer type output with header
        IPD   = 1
        IF (NREF.EQ.PRINTF) CALL TABHED ('SWASH', PRINTF)                 30.20
      ELSE IF (RTYPE .EQ. 'BLKD') THEN
!       output to datafile without header
        IPD = 2
      ELSE
        IPD = 3
      ENDIF
      IF (ITEST.GE.90) WRITE (PRTEST, 21)  RTYPE,NREF, OQI(3)             40.31 40.03
  21  FORMAT (' Test SWBLOK: RTYPE NREF NVAR ',A4,2(1X,I6))
      FILENM = OUTP_FILES(OQI(2))                                         40.31 40.13
      MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.                          40.41 40.30
     &         INDEX (FILENM, '.mat' ).NE.0                               40.41 40.30
      IF (NREF.EQ.0) THEN
        IOSTAT = -1                                                       20.75
        CALL FOR (NREF, FILENM, 'UF', IOSTAT)
        IF (STPNOW()) RETURN                                              34.01
        OQI(1) = NREF                                                     40.31 30.00
        OUTP_FILES(OQI(2)) = FILENM                                       40.41
        IF (MATLAB) THEN                                                  40.30
           CLOSE(NREF)                                                    40.30
           OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',               40.30
     &          STATUS='REPLACE',
!MatL4     &          ACCESS='DIRECT', RECL=1)                                  40.30
     &          ACCESS='DIRECT', RECL=4)                                  41.08 40.30
           IREC(IRQ) = 1                                                  40.51
        END IF
      ENDIF
      IDLA = OQI(4)                                                       40.31 30.00
      NVAR = OQI(3)                                                       40.31 30.00
!
      IF (ITEST.GE.90) WRITE (PRTEST, 22)  NREF, FILENM
  22  FORMAT (' Test SWBLOK: NREF FILENM  ', I6, A40)
!
      CTIM = CHTIME                                                       40.41
      CALL TXPBLA(CTIM,IF,IL)                                             40.41
      IF (ITMOPT.EQ.1) THEN
         CTIM(9:9)='_'                                                    40.41
      ELSE IF (ITMOPT.EQ.7) THEN
         CTIM(7:7)='_'
      ENDIF
!
      IF ( KMAX.LT.10 ) THEN
         CRFORM = '(I1.1)'
         KFML   = 56
      ELSE IF ( KMAX.LT.100 ) THEN
         CRFORM = '(I2.2)'
         KFML   = 55
      ELSE IF ( KMAX.LT.1000 ) THEN
         CRFORM = '(I3.3)'
         KFML   = 54
      ELSE
         CRFORM = '(I4.4)'
         KFML   = 53
      ENDIF
!
      DO  800  JVAR=1,NVAR
        IVTYPE = IVTYP(JVAR)                                              40.31 30.00
        DFAC   = FAC(JVAR)                                                40.31 30.00
!
        IF ( IVTYPE.LT. 4 .OR. IVTYPE.EQ. 5 .OR. IVTYPE.GT.20 .AND.
     &       IVTYPE.LT.27 .OR. IVTYPE.GT.32 .AND. IVTYPE.LT.38 ) THEN
           IVSTAT1 = .TRUE.
        ELSE
           IVSTAT1 = .FALSE.
        ENDIF
        IF ( IVTYPE.GT.56 .AND. IVTYPE.LT.60 ) THEN
           IVSTAT2 = .TRUE.
        ELSE
           IVSTAT2 = .FALSE.
        ENDIF
        IF ( IVTYPE.GT.83 .AND. IVTYPE.NE.89 .AND.
     &       IVTYPE.NE.90 .AND. IVTYPE.NE.91      ) THEN
           IVSTAT3 = .TRUE.
        ELSE
           IVSTAT3 = .FALSE.
        ENDIF
!
        IF (IPD.EQ.1) THEN
          IF (DFAC.LE.0.) THEN
!           determine default factor for print output
            IF (OVHEXP(IVTYPE) .LT. 0.5E10) THEN
              IFAC = INT (10.+LOG10(OVHEXP(IVTYPE))) - 13                 30.20
            ELSE
              FMAX = 1000.
              IF ( IVTYPE.LE.50 ) THEN
                 IF (OVSVTY(IVTYPE).EQ.1) THEN
                    FMAX = MAXVAL(ABS(VOQ(:,VOQR(IVTYPE))))
                 ELSE IF (OVSVTY(IVTYPE).EQ.3) THEN
                    FTIP1 = MAXVAL(ABS(VOQ(:,VOQR(IVTYPE))))
                    FTIP2 = MAXVAL(ABS(VOQ(:,VOQR(IVTYPE)+1)))
                    FMAX  = MAX(FTIP1, FTIP2)
                 ENDIF
              ELSE
                 IF (OVSVTY(IVTYPE).EQ.1) THEN
                    FMAX = MAXVAL(ABS(VOQK(:,:,VOQR(IVTYPE))))
                 ELSE IF (OVSVTY(IVTYPE).EQ.3) THEN
                    FTIP1 = MAXVAL(ABS(VOQK(:,:,VOQR(IVTYPE))))
                    FTIP2 = MAXVAL(ABS(VOQK(:,:,VOQR(IVTYPE)+1)))
                    FMAX  = MAX(FTIP1, FTIP2)
                 ENDIF
              ENDIF
              IFAC = INT (10.+LOG10(FMAX)) - 13
            ENDIF
            DFAC = 10.**IFAC
          ENDIF
        ELSE
          IF (DFAC.LE.0.) DFAC = 1.
        ENDIF
!
        IF (ITEST .GE. 80) WRITE(PRTEST, 6020) JVAR, IVTYPE, DFAC,
     &    COSCQ, SINCQ
 6020   FORMAT(' Test SWBLOK: jvar, ivtype, dfac, coscq, sincq',
     &          2I10,3E12.5)
!
        IF (OVSVTY(IVTYPE) .LT. 3) THEN
!                      scalar quantities
          IF (MATLAB) THEN                                                40.30
             IF ( IVSTAT1 ) THEN                                          40.94 40.41
                NAMVAR = OVSNAM(IVTYPE)                                   40.41
                CALL SWRMAT( MYK, MXK, NAMVAR,                            40.41
     &                       VOQ(1,VOQR(IVTYPE)), NREF, IREC(IRQ),        40.51
     &                       IDLA, OVEXCV(IVTYPE) )
             ELSE IF ( IVSTAT2 ) THEN
                DO K = 0, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      KNUM
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                          VOQK(1,K,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                          IDLA, OVEXCV(IVTYPE) )
                ENDDO
             ELSE IF ( IVSTAT3 ) THEN
                DO K = 1, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      KNUM
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                          VOQK(1,K,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                          IDLA, OVEXCV(IVTYPE) )
                ENDDO
             ELSE IF ( IVTYPE.LE.50 ) THEN                                40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))       40.41
                IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                '_'//CTIM
                CALL SWRMAT( MYK, MXK, NAMVAR,                            40.41
     &                       VOQ(1,VOQR(IVTYPE)), NREF, IREC(IRQ),        40.51
     &                       IDLA, OVEXCV(IVTYPE) )
             ELSE IF ( IVTYPE.LE.70 ) THEN
                DO K = 0, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      KNUM
                   IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                   '_'//CTIM
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                          VOQK(1,K,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                          IDLA, OVEXCV(IVTYPE) )
                ENDDO
             ELSE
                DO K = 1, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      KNUM
                   IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                   '_'//CTIM
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                          VOQK(1,K,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                          IDLA, OVEXCV(IVTYPE) )
                ENDDO
             END IF                                                       40.41
          ELSE
             IF ( IVTYPE.LE.50 ) THEN
                CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &          MXK, MYK, IDLA, OVLNAM(IVTYPE), VOQ(1,VOQR(IVTYPE)))
             ELSE IF ( IVTYPE.LE.70 ) THEN
                DO K = 0, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                        MXK, MYK, IDLA,
     &                        OVLNAM(IVTYPE)(1:KFML)//' k='//TRIM(KNUM),
     &                        VOQK(1,K,VOQR(IVTYPE)))
                ENDDO
             ELSE
                DO K = 1, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                        MXK, MYK, IDLA,
     &                        OVLNAM(IVTYPE)(1:KFML)//' k='//TRIM(KNUM),
     &                        VOQK(1,K,VOQR(IVTYPE)))
                ENDDO
             ENDIF
          END IF
        ELSE
!                     vectorial quantities
          IF (MATLAB) THEN                                                40.30
             IF ( IVSTAT1 ) THEN
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                   '_x'
                CALL SWRMAT( MYK, MXK, NAMVAR,
     &                    VOQ(1,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                    IDLA, OVEXCV(IVTYPE))
             ELSE IF ( IVSTAT2 ) THEN
                DO K = 0, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_x'
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                       VOQK(1,K,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                       IDLA, OVEXCV(IVTYPE))
                ENDDO
             ELSE IF ( IVSTAT3 ) THEN
                DO K = 1, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_x'
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                       VOQK(1,K,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                       IDLA, OVEXCV(IVTYPE))
                ENDDO
             ELSE IF ( IVTYPE.LE.50 ) THEN                                40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//     40.41
     &                   '_x'
                IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                '_'//CTIM
                CALL SWRMAT( MYK, MXK, NAMVAR,                            40.41
     &                    VOQ(1,VOQR(IVTYPE)), NREF, IREC(IRQ),           40.51
     &                    IDLA, OVEXCV(IVTYPE))
             ELSE IF ( IVTYPE.LE.70 ) THEN
                DO K = 0, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_x'
                   IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                   '_'//CTIM
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                       VOQK(1,K,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                       IDLA, OVEXCV(IVTYPE))
                ENDDO
             ELSE
                DO K = 1, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_x'
                   IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                   '_'//CTIM
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                       VOQK(1,K,VOQR(IVTYPE)), NREF, IREC(IRQ),
     &                       IDLA, OVEXCV(IVTYPE))
                ENDDO
             END IF                                                       40.41
             IF ( IVSTAT1 ) THEN                                          40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//     40.41
     &                   '_y'
                CALL SWRMAT( MYK, MXK, NAMVAR,                            40.41
     &                    VOQ(1,VOQR(IVTYPE)+1), NREF, IREC(IRQ),         40.51
     &                    IDLA, OVEXCV(IVTYPE))
             ELSE IF ( IVSTAT2 ) THEN
                DO K = 0, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_y'
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                       VOQK(1,K,VOQR(IVTYPE)+1), NREF, IREC(IRQ),
     &                       IDLA, OVEXCV(IVTYPE))
                ENDDO
             ELSE IF ( IVSTAT3 ) THEN
                DO K = 1, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_y'
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                       VOQK(1,K,VOQR(IVTYPE)+1), NREF, IREC(IRQ),
     &                       IDLA, OVEXCV(IVTYPE))
                ENDDO
             ELSE IF ( IVTYPE.LE.50 ) THEN                                40.41
                NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//     40.41
     &                   '_y'
                IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                '_'//CTIM
                CALL SWRMAT( MYK, MXK, NAMVAR,                            40.41
     &                    VOQ(1,VOQR(IVTYPE)+1), NREF, IREC(IRQ),         40.51
     &                    IDLA, OVEXCV(IVTYPE))
             ELSE IF ( IVTYPE.LE.70 ) THEN
                DO K = 0, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_y'
                   IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                   '_'//CTIM
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                       VOQK(1,K,VOQR(IVTYPE)+1), NREF, IREC(IRQ),
     &                       IDLA, OVEXCV(IVTYPE))
                ENDDO
             ELSE
                DO K = 1, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_y'
                   IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                   '_'//CTIM
                   CALL SWRMAT( MYK, MXK, NAMVAR,
     &                       VOQK(1,K,VOQR(IVTYPE)+1), NREF, IREC(IRQ),
     &                       IDLA, OVEXCV(IVTYPE))
                ENDDO
             END IF                                                       40.41
          ELSE
             IF ( IVTYPE.LE.50 ) THEN
                CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &          MXK, MYK, IDLA,
     &          OVLNAM(IVTYPE)(1:52)//'  X-comp', VOQ(1,VOQR(IVTYPE)))
                CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &          MXK, MYK, IDLA,
     &          OVLNAM(IVTYPE)(1:52)//'  Y-comp', VOQ(1,VOQR(IVTYPE)+1))
             ELSE IF ( IVTYPE.LE.70 ) THEN
                DO K = 0, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &             MXK, MYK, IDLA,
     &             OVLNAM(IVTYPE)(1:KFML-8)//' k='//TRIM(KNUM)//
     &             '  X-comp',
     &             VOQK(1,K,VOQR(IVTYPE)))
                   CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &             MXK, MYK, IDLA,
     &             OVLNAM(IVTYPE)(1:KFML-8)//' k='//TRIM(KNUM)//
     &             '  Y-comp',
     &             VOQK(1,K,VOQR(IVTYPE)+1))
                ENDDO
             ELSE
                DO K = 1, KMAX
                   WRITE(KNUM, FMT=CRFORM) K
                   CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &             MXK, MYK, IDLA,
     &             OVLNAM(IVTYPE)(1:KFML-8)//' k='//TRIM(KNUM)//
     &             '  X-comp',
     &             VOQK(1,K,VOQR(IVTYPE)))
                   CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &             MXK, MYK, IDLA,
     &             OVLNAM(IVTYPE)(1:KFML-8)//' k='//TRIM(KNUM)//
     &             '  Y-comp',
     &             VOQK(1,K,VOQR(IVTYPE)+1))
                ENDDO
             ENDIF
          END IF
        ENDIF
!
  800 CONTINUE
      IF (IPD.EQ.1 .AND. NREF.EQ.PRINTF) WRITE (PRINTF, 6030)
 6030 FORMAT (///)
!
      RETURN
! * end of subroutine SWBLOK *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SBLKPT (IPD, NREF, DFAC, PSNAME, QUNIT,
     &                    MXK, MYK, IDLA, STRING, OQVALS)
!                                                                      *
!***********************************************************************

      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE SwashCommdata3                                                  40.41
      USE outp_data                                                       40.13
      USE SwashTimeComm                                                   40.41
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWASH team                               |
!   --|-----------------------------------------------------------|--
!
!
!     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
!     Copyright (C) 2010-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.82: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     00.00, Mar. 87: subroutine heading added, some variable names
!                     line numbers changed, layout modified
!     00.04, Feb. 90: lay-out of output changed according to IDLA=1
!     30.72, Sept 97: Changed DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.82, Nov. 98: Corrected syntax format statement
!     40.13, July 01: variable formats introduced, using module OUTP_DATA
!     40.13, Oct. 01: longer output filenames now obtained from array
!                     OUTP_FILES (in module OUTP_DATA)
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Writing the block output either on paper or to datafile
!
!  3. Method
!
!     ---
!
!  4. PARAMETERLIST
!
!     IPD     INT    input    switch for printing on paper (IPD=1)
!                             or writing to datafile (IPD = 2 or 3)
!     NREF    INT    input    unit reference number of output file
!     DFAC    REAL   input    multiplication factor of block output
!     IVTYPE  INT    input    type of the output quantity
!                             Note: IVTYPE=0 for Y-component of a
!                             vectorial quantity
!     PSNAME  CH*8   input    name of output point set (frame)
!     QUNIT   CH*6   input    physical unit (dimension) of variable
!     MXK     int    input    number of points in x-direction of frame
!     MYK     int    input    number of points in y-direction of frame
!     IDLA    INT    input    controls lay-out of output (see user manual)
!     STRING  CH*(*) input    description of output variable
!
!  8. Subroutines used
!
!       ---
!
!  9. Subroutines calling
!
!       SWBLOK (SWAN/OUTP)
!
! 10. Error messages
!
!       ---
!
! 11. Remarks
!
!       ---
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       If IPD = 1 (output on paper) then
!           If DFAC < 0 (DFAC not given by the user) then
!               Compute maximum value of output variable
!               Compute multiplication factor DFAC
!           ------------------------------------------------------------
!           Print block heading
!           For each IX of the output frame do
!               Print IX and for every IY the value of the outputvariable
!           ------------------------------------------------------------
!       Else
!           If DFAC < 0 then DFAC = 1.
!           Write output variable line by line to datafile
!       ----------------------------------------------------------------
!
! 13. Source text
!

      CHARACTER (LEN=20) :: WFORM1 = '(A1, 2X, 151(I6))'                  40.13
      CHARACTER (LEN=21) :: WFORM2 = '(1X,I3,1X, 151(F6.0))'              40.13
      CHARACTER (LEN=20) :: WFORM3 = '(5X, 151(F6.0))'                    40.13

      CHARACTER PSNAME*8, STRING*(*), QUNIT*(*)                           40.00
      REAL      DFAC, OQVALS(*)
      INTEGER   NREF, MXK, MYK, IPD
      LOGICAL   BPRN
      SAVE IENT
      DATA IENT /0/
      IF (LTRACE) CALL STRACE (IENT,'SBLKPT')
!
      IF (ITEST.GE.150) WRITE (PRTEST, 10) NREF,IPD,MXK,MYK
  10  FORMAT (' SBLKPT', 4(I6))                                           30.82
!
!
!     divide all output values by the given factor (DFAC)
!
      IF (ABS(DFAC-1.) .GT. 0.001) THEN
        RPDFAC=1./DFAC
        DO 15 IP = 1, MXK*MYK
           OQVALS(IP) = OQVALS(IP)*RPDFAC
   15   CONTINUE
      ENDIF
!
!
!      IFF = VOQR(IVTYPE)
!
      IF (IPD.EQ.1) THEN
!
!       ***** output on paper *****
!
        WRITE (NREF, 20) OUT_COMMENT                                      40.13
        WRITE (NREF, 20) OUT_COMMENT                                      40.13
  20    FORMAT (A)
        WRITE (NREF, 22) OUT_COMMENT, PROJNR, PSNAME, STRING,             40.13
     &                   DFAC, QUNIT                                      40.13
  22    FORMAT (A,' Run:', A4, '  Frame:  ',A8,' **  ',A,', Unit:',       40.13
     &          E12.4, 1X, A)                                             40.13
        IF (NSTATM .GT. 0) THEN
          WRITE (NREF, 24) OUT_COMMENT, CHTIME                            40.13
  24      FORMAT (A,' Time:', A)                                          40.13
        ELSE
          WRITE (NREF, 20) OUT_COMMENT                                    40.13
        ENDIF                                                             40.13
        WRITE (NREF, 20) OUT_COMMENT                                      40.13

        ISP = 151                                                         30.21
        DO  31  IXP1 = 1, MXK, ISP                                        30.72
          IXP2 = IXP1+ISP-1
          IF (IXP2.GT.MXK) IXP2=MXK

          WRITE (WFORM1(15:15), '(I1)') DEC_BLOCK                         40.13
          WRITE (WFORM2(17:17), '(I1)') DEC_BLOCK                         40.13
          WRITE (WFORM3(11:11), '(I1)') DEC_BLOCK                         40.13
          IF (ITEST.GE.80) WRITE (PRTEST, 25) WFORM1, WFORM2, WFORM3      40.13
  25      FORMAT (' SBLKPT Formats: ', A, /, 17X, A, /, 17X, A)           40.13

          WRITE (NREF, 26) OUT_COMMENT                                    40.13
  26      FORMAT (A1,'         X --->')                                   40.13
          WRITE (NREF, 20) OUT_COMMENT                                    40.13
          WRITE (NREF, WFORM1) OUT_COMMENT, (II-1,II=IXP1,IXP2)           40.13
          WRITE (NREF, 99030) OUT_COMMENT                                 40.13
99030     FORMAT (A1, 'Y')                                                40.13

          BPRN = .TRUE.
          DO 30 IYK = MYK, 1, -1
            IP = (IYK-1)*MXK
            IF (BPRN) THEN
              WRITE (NREF, WFORM2) IYK-1,
     &        (OQVALS(IP+IXK), IXK=IXP1,IXP2)
            ELSE
              WRITE (NREF, WFORM3)
     &        (OQVALS(IP+IXK), IXK=IXP1,IXP2)
            ENDIF
!!!            BPRN = .NOT. BPRN
   30     CONTINUE                                                        30.72
   31   CONTINUE                                                          30.72
      ELSE
!
!       ***** output to datafile *****
!
        ISP=6
        IF (IDLA.EQ.4) THEN
          WRITE (NREF, FLT_BLOCK) (OQVALS(IP), IP=1, MXK*MYK)             40.13
        ELSE
          DO 50 IYK = 1, MYK                                              13/FEB
            IF (IDLA.EQ.3) THEN
              IP = (IYK-1)*MXK
            ELSE
              IP = (MYK-IYK)*MXK
            ENDIF
            WRITE (NREF, FLT_BLOCK) (OQVALS(IP+IXK), IXK=1,MXK)           40.13
   50     CONTINUE
        ENDIF
      ENDIF
!
      RETURN
! * end of subroutine SBLKPT *
      END
!****************************************************************
!
      SUBROUTINE SWBLKP ( OQI   , IVTYP, MXK  , MYK, VOQR, VOQ,
     &                    VOQK  , IONOD )                                 40.51
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE SwashCommdata3, ONLY: KMAX
      USE outp_data
      USE M_PARALL                                                        40.51
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
!     Copyright (C) 2010-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!
!  1. Updates
!
!     40.31, Dec. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: further optimization
!
!  2. Purpose
!
!     Write block data to process output file
!
!  3. Method
!
!     Write data without header using the format FLT_BLKP
!
!  4. Argument variables
!
!     IONOD       array indicating in which subdomain output points       40.51
!                 are located                                             40.51
!     IVTYP       type of variable output
!     MXK         number of points in x-direction of output frame
!     MYK         number of points in y-direction of output frame
!     OQI         array containing output request data
!     VOQ         output variables
!     VOQR        array containing information for output
!
      INTEGER MXK, MYK, OQI(4), IVTYP(OQI(3)), VOQR(*)
      INTEGER IONOD(*)                                                    40.51
      REAL    VOQ(MXK*MYK,*), VOQK(MXK*MYK,0:KMAX,*)
!
!  6. Local variables
!
!     IENT  :     number of entries
!     IOSTAT:     status of input/output
!     IP    :     pointer
!     IVTYPE:     type number output variable
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     K     :     counter
!     NREF  :     unit reference number
!     NVAR  :     number of variables
!
      INTEGER IENT, IOSTAT, IP, IVTYPE, IXK, IYK, JVAR, K, NREF, NVAR
!
!  8. Subroutines used
!
!     FOR              General open file routine
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWOUTP
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBLKP')

      NREF = OQI(1)
      IF (NREF.EQ.0) THEN
         FILENM = OUTP_FILES(OQI(2))
         IOSTAT = -1
         CALL FOR (NREF, FILENM, 'UU', IOSTAT)
         IF (STPNOW()) RETURN
         OQI(1) = NREF
         OUTP_FILES(OQI(2)) = FILENM                                      40.41
      END IF
      NVAR = OQI(3)

      DO JVAR = 1, NVAR
         IVTYPE = IVTYP(JVAR)
         IF ( IVTYPE.LE.50 ) THEN
            DO IYK = 1, MYK                                               40.51
               IP = (IYK-1)*MXK                                           40.51
               DO IXK = 1, MXK                                            40.51
                  IF ( IONOD(IP+IXK).EQ.INODE )                           40.51
     &               WRITE (NREF) VOQ(IP+IXK,VOQR(IVTYPE))                40.51
               END DO                                                     40.51
            END DO                                                        40.51
            IF ( OVSVTY(IVTYPE).GE.3 ) THEN
               DO IYK = 1, MYK                                            40.51
                  IP = (IYK-1)*MXK                                        40.51
                  DO IXK = 1, MXK                                         40.51
                     IF ( IONOD(IP+IXK).EQ.INODE )                        40.51
     &                  WRITE (NREF) VOQ(IP+IXK,VOQR(IVTYPE)+1)           40.51
                  END DO                                                  40.51
               END DO                                                     40.51
            END IF
         ELSEIF ( IVTYPE.LE.70 ) THEN
            DO IYK = 1, MYK
               IP = (IYK-1)*MXK
               DO IXK = 1, MXK
                  IF ( IONOD(IP+IXK).EQ.INODE )
     &               WRITE (NREF)
     &                           (VOQK(IP+IXK,K,VOQR(IVTYPE)), K=0,KMAX)
               END DO
            END DO
            IF ( OVSVTY(IVTYPE).GE.3 ) THEN
               DO IYK = 1, MYK
                  IP = (IYK-1)*MXK
                  DO IXK = 1, MXK
                     IF ( IONOD(IP+IXK).EQ.INODE )
     &                  WRITE (NREF)
     &                         (VOQK(IP+IXK,K,VOQR(IVTYPE)+1), K=0,KMAX)
                  END DO
               END DO
            END IF
         ELSE
            DO IYK = 1, MYK
               IP = (IYK-1)*MXK
               DO IXK = 1, MXK
                  IF ( IONOD(IP+IXK).EQ.INODE )
     &               WRITE (NREF)
     &                           (VOQK(IP+IXK,K,VOQR(IVTYPE)), K=1,KMAX)
               END DO
            END DO
            IF ( OVSVTY(IVTYPE).GE.3 ) THEN
               DO IYK = 1, MYK
                  IP = (IYK-1)*MXK
                  DO IXK = 1, MXK
                     IF ( IONOD(IP+IXK).EQ.INODE )
     &                  WRITE (NREF)
     &                         (VOQK(IP+IXK,K,VOQR(IVTYPE)+1), K=1,KMAX)
                  END DO
               END DO
            END IF
         ENDIF
      END DO

      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWTABP (RTYPE , OQI  , IVTYP, PSNAME, MIP, VOQR, VOQ,    40.31
     &                   VOQK  , IONOD)                                   40.51
!                                                                      *
!***********************************************************************

      USE OCPCOMM2                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE SwashCommdata2                                                  40.41
      USE SwashCommdata3                                                  40.41
      USE SwashCommdata4                                                  40.41
      USE outp_data                                                       40.13
      USE SwashTimeComm                                                   40.41
      USE M_PARALL                                                        40.51
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWASH team                               |
!   --|-----------------------------------------------------------|--
!
!
!     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
!     Copyright (C) 2010-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.62: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.80: Nico Booij
!     30.82: IJsbrand Haagsma
!     32.01: Roeland Ris & Cor van der Schelde
!     34.01: Jeroen Adema
!     40.00: Nico Booij (Non-stationary boundary conditions)
!     40.03, 40.13: Nico Booij
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.51: Agnieszka Herman
!     SWASH: Marcel Zijlema
!
!  1. Updates
!
!     30.50, Sep. 96: option TABI (indexed file) added
!     30.62, Jul. 97: corrected initialisation of table output
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     32.01, Jan. 98: Extended initialisation of NUMDEC for SETUP
!     30.80, Apr. 98: number of decimals for setup from 2 to 3
!     40.00, June 98: severely revised
!     30.82, Oct. 98: Header information is now also printed in PRINT file
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Mar. 00: number of decimals (NUMDEC) is made larger
!     40.13, Jan. 01: program version now written into table heading
!            Mar. 01: XOFFS and YOFFS were incorrectly added to coordinates
!                     (they are already included in VOQ values)
!     40.13, July 01: variable formats introduced, using module OUTP_DATA
!                     comment sign in front of heading lines
!     40.13, Oct. 01: longer output filenames now obtained from array
!                     OUTP_FILES (in module OUTP_DATA)
!     40.31, Dec. 03: removing POOL construction
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: further optimization
!     SWASH, Sep. 10: extension with layer-dependent quantities
!
!  2. Purpose
!
!     Printing of output in the form of a table for any type of
!     output point set
!
!  3. Method
!
!     A table is made in which for each point the required output
!     variables are printed in the order given by the user. If more
!     variables are required than one line can contain, writing is
!     continued on the next line before output for the next point is
!     started.
!
!  4. Argument variables
!
!     PSNAME
! i   RTYPE : Type of output request
!             ='TABD'; Output to datafile (no header information)
!             ='TABI'; Indexed output for table in ArcView format
!             ='TABP'; Output to paper (with header information)
!             ='TABS';
!             ='TABT';
!
      CHARACTER RTYPE*4, PSNAME*8
!
!     MIP
!     VOQR
!
      INTEGER   MIP, VOQR(*), OQI(4), IVTYP(OQI(3))
      INTEGER   IONOD(*)                                                  40.51
!
!     VOQ
!
      REAL      VOQ(MIP,*), VOQK(MIP,0:KMAX,*)
!
!  5. Parameter variables
!
!     MXOUTL
!
      INTEGER MXOUTL
!
      PARAMETER (MXOUTL=720)
!
!  6. Local variables
!
!     NUMDEC
!
      INTEGER NUMDEC
!
!  8. Subroutines used
!
!     SUHEAD (all SWAN/OUTP)
!     FOR
!     TABHED (all Ocean Pack)
!
      LOGICAL STPNOW                                                      34.01
!
!  9. Subroutines calling
!
!     OUTPUT (SWAN/OUTP)
!
! 10. ERROR MESSAGES
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       If unit ref. number = 0
!       Then Read filename from array IOUTR
!            Call FOR to open datafile
!            If Rtype = 'TABP' or 'TABI'
!            Then Print heading for required table
!       ----------------------------------------------------------------
!       If Rtype = 'TABS'
!       Then write time into file
!       ----------------------------------------------------------------
!       Make Output line blank
!       Make Linkar = 1
!       If Rtype = 'TABI'
!       Then write index into output line
!            update Linkar
!       ----------------------------------------------------------------
!       For every output point do
!           For all output quantities do
!               Get value for array VOQ
!               If output quantity is TIME
!               Then make Format='(A18)'
!                    write time into output line
!                    make lfield = 18
!               Else if Rtype = 'TABD'
!                    Then Format = '(E12.4)'
!                         make lfield = 12
!                    Else Format = '(F13.X)'
!                         make lfield = 13
!                         determine number of decimals and write into
!                         Format
!                    ---------------------------------------------------
!                    Write value into output line according to Format
!               --------------------------------------------------------
!               Make Linkar = Linkar + lfield + 1
!           ------------------------------------------------------------
!           Write Output line to file
!       ----------------------------------------------------------------
!
! 13. Source text
!
      CHARACTER (LEN=15)     :: FSTR                                      40.22 40.13
      CHARACTER (LEN=MXOUTL) :: OUTLIN                                    40.13
      CHARACTER (LEN=8)      :: CRFORM = '(2F14.4)'                       40.03

      INTEGER, SAVE :: IENT=0                                             40.13
      IF (LTRACE) CALL STRACE(IENT,'SWTABP')
!
      NREF = OQI(1)                                                       40.31 30.00
      NVAR = OQI(3)                                                       40.31 30.00
!
!     Header information is printed once for each data file and for
!     each entry in this routine when the table is written to the PRINT file
!
      IF (NREF .EQ. 0 .OR. NREF.EQ.PRINTF) THEN                           30.82
        IF (NREF.EQ.0) THEN                                               30.82
          FILENM = OUTP_FILES(OQI(2))                                     40.31 40.13
          IOSTAT = -1                                                     20.75
          CALL FOR (NREF, FILENM, 'UF', IOSTAT)
          IF (STPNOW()) RETURN                                            34.01
          OQI(1) = NREF                                                   40.31 30.00
          OUTP_FILES(OQI(2)) = FILENM                                     40.41
        END IF                                                            30.82
        IF (RTYPE .NE. 'TABD') THEN
          OUTLIN = '    '
!
!         write heading into file
!
          IF (RTYPE.EQ.'TABP' .OR. RTYPE.EQ.'TABI') THEN
            WRITE (NREF, 20) OUT_COMMENT                                  40.13
            WRITE (NREF, 20) OUT_COMMENT                                  40.13
  20        FORMAT (A)
            WRITE (NREF, 43) OUT_COMMENT, PROJNR, PSNAME, TRIM(VERTXT)    40.13
  43        FORMAT (A1, ' Run:', A4,'  Table:',A8, 10X,                   40.13
     &      'SWASH version:', A)                                          40.13
            WRITE (NREF, 20) OUT_COMMENT                                  40.13
!           write (short) names of output quantities
            IF (RTYPE.EQ.'TABI') THEN
              OUTLIN(1:12) = OUT_COMMENT // '           '                 40.13
              LINKAR = 12
            ELSE
              OUTLIN(1:3) = OUT_COMMENT                                   40.13
              LINKAR = 4                                                  40.13
            ENDIF
            DO  JVAR = 1, NVAR
              IVTYPE = IVTYP(JVAR)                                        40.31
              IF (IVTYPE.EQ.40) THEN
                LFIELD = 18
              ELSE
                LFIELD = 13
              ENDIF
              IF (OVSVTY(IVTYPE).LE.2) THEN
                OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '     '//OVSNAM(IVTYPE)//'              '
              ELSE
                OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '   X-'//OVSNAM(IVTYPE)//'              '
                LINKAR = LINKAR+LFIELD+1
                OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '   Y-'//OVSNAM(IVTYPE)//'              '
              ENDIF
              LINKAR = LINKAR+LFIELD+1
            ENDDO
            WRITE (NREF, '(A)') OUTLIN(1:LINKAR-1)
!           write units of output quantities
            OUTLIN = '    '
            IF (RTYPE.EQ.'TABI') THEN
              OUTLIN(1:12) = OUT_COMMENT // '           '                 40.13
              LINKAR = 12
            ELSE
              OUTLIN(1:3) = OUT_COMMENT                                   40.13
              LINKAR = 4                                                  40.13
            ENDIF
            DO  JVAR = 1, NVAR
              IVTYPE = IVTYP(JVAR)                                        40.31
              IF (IVTYPE.EQ.40) THEN
                LFIELD = 18
              ELSE
                LFIELD = 13
              ENDIF
              DO ISTR = LEN(OVUNIT(IVTYPE)), 1, -1
                IF (OVUNIT(IVTYPE)(ISTR:ISTR) .NE. ' ') THEN
                  LSTR = ISTR
                  GOTO 51
                ENDIF
              ENDDO
              LSTR = 1
  51          OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '     ['//OVUNIT(IVTYPE)(1:LSTR)//']            '
              IF (OVSVTY(IVTYPE).GT.2) THEN
                LINKAR = LINKAR+LFIELD+1
                OUTLIN(LINKAR:LINKAR+LFIELD) =
     &                '     ['//OVUNIT(IVTYPE)(1:LSTR)//']            '
              ENDIF
              LINKAR = LINKAR+LFIELD+1
            ENDDO
            WRITE (NREF, '(A)') OUTLIN(1:LINKAR-1)
            WRITE (NREF, 20) OUT_COMMENT                                  40.13
          ELSE IF (RTYPE.EQ.'TABT' .OR. RTYPE.EQ.'TABS' .OR.
     &             RTYPE.EQ.'TABV' .OR. RTYPE.EQ.'TABU') THEN
            WRITE (NREF, 101) 1
 101        FORMAT ('SWASH', I4, T41, 'SWASH standard file, version')
            WRITE (NREF, 111) OUT_COMMENT, VERTXT                         40.13
 111        FORMAT (A1, '   Data produced by SWASH version ', A)          40.13
            WRITE (NREF, 113) OUT_COMMENT, PROJID, PROJNR                 40.13
 113        FORMAT (A1, '   Project: ', A, ';  run number: ', A)          40.13
 102        FORMAT (A, T41, A)
 103        FORMAT (I6, T41, A)
            IF (RTYPE.EQ.'TABT' .OR. RTYPE.EQ.'TABV') THEN
              WRITE (NREF,102) 'TABLE'
            ELSE
              IF (NSTATM.EQ.1) THEN
                WRITE (NREF, 102) 'TIME', 'time-dependent data'
                WRITE (NREF, 103) ITMOPT, 'time coding option'
              ENDIF
              IF (KSPHER.EQ.0) THEN                                       33.09
                WRITE (NREF, 102) 'LOCATIONS', 'locations in x-y-space'
                CRFORM = '(2F14.4)'                                       40.03
              ELSE                                                        33.09
                WRITE (NREF, 102) 'LONLAT',
     &                    'locations in spherical coordinates'            33.09
                CRFORM = '(2F12.6)'                                       40.13
              ENDIF                                                       33.09
              WRITE (NREF, 103) MIP, 'number of locations'
              DO 110 IP = 1, MIP
                WRITE (NREF, FMT=CRFORM)                                  40.03
     &                DBLE(VOQ(IP,VOQR(1))), DBLE(VOQ(IP,VOQR(2)))        40.13
 110          CONTINUE
            ENDIF
            WRITE (NREF, 102) 'QUANT', 'description of quantities'
            NKOLS = NVAR
            DO  JVAR = 1, NVAR
              IVTYPE = IVTYP(JVAR)                                        40.31
              IF (OVSVTY(IVTYPE).GT.2) NKOLS = NKOLS + 1
            ENDDO
            WRITE (NREF, 103) NKOLS, 'number of quantities in table'
            DO  JVAR = 1, NVAR
              IVTYPE = IVTYP(JVAR)                                        40.31
              IF (OVSVTY(IVTYPE).LE.2) THEN
                WRITE (NREF, 102) OVSNAM(IVTYPE), OVLNAM(IVTYPE)
                WRITE (NREF, 102) OVUNIT(IVTYPE), 'unit'
                IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 ) THEN
                   WRITE (NREF, 103) 1, 'dimension in vertical'
                ELSE IF ( IVTYPE.LE.70 ) THEN
                   WRITE (NREF, 103) KMAX+1, 'dimension in vertical'
                ELSE
                   WRITE (NREF, 103) KMAX, 'dimension in vertical'
                ENDIF
                WRITE (NREF, 104) OVEXCV(IVTYPE), 'exception value'
 104            FORMAT (E14.4, T41, A)
              ELSE
                WRITE (NREF, 102) 'X-'//OVSNAM(IVTYPE), OVLNAM(IVTYPE)
                WRITE (NREF, 102) OVUNIT(IVTYPE), 'unit'
                IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 ) THEN
                   WRITE (NREF, 103) 1, 'dimension in vertical'
                ELSE IF ( IVTYPE.LE.70 ) THEN
                   WRITE (NREF, 103) KMAX+1, 'dimension in vertical'
                ELSE
                   WRITE (NREF, 103) KMAX, 'dimension in vertical'
                ENDIF
                WRITE (NREF, 104) OVEXCV(IVTYPE), 'exception value'
                WRITE (NREF, 102) 'Y-'//OVSNAM(IVTYPE), OVLNAM(IVTYPE)
                WRITE (NREF, 102) OVUNIT(IVTYPE), 'unit'
                IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 ) THEN
                   WRITE (NREF, 103) 1, 'dimension in vertical'
                ELSE IF ( IVTYPE.LE.70 ) THEN
                   WRITE (NREF, 103) KMAX+1, 'dimension in vertical'
                ELSE
                   WRITE (NREF, 103) KMAX, 'dimension in vertical'
                ENDIF
                WRITE (NREF, 104) OVEXCV(IVTYPE), 'exception value'
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
!
!     ***** printing of the table *****
!
      IF (RTYPE.EQ.'TABS' .OR. RTYPE.EQ.'TABU') THEN
        IF (ITMOPT.NE.7) THEN
           OUTLIN = 'date and time'
        ELSE
           OUTLIN = 'time'
        ENDIF
        IF (NSTATM.EQ.1) WRITE (NREF, 102) CHTIME, TRIM(OUTLIN)
      ENDIF
      DO 70 IP = 1, MIP
       IF ( .NOT.PARLL .OR. IONOD(IP).EQ.INODE ) THEN                     40.51
        LINKAR = 1
        OUTLIN = '    '
        IF (RTYPE.EQ.'TABI') THEN
!         write point sequence number as first column
          WRITE (OUTLIN(1:8), '(I8)') IP
          LINKAR = 9
          OUTLIN(LINKAR:LINKAR) = ' '
        ENDIF
        DO 60 JVAR = 1, NVAR
          IVTYPE = IVTYP(JVAR)                                            40.31
          IF (IVTYPE.EQ.40) THEN                                          40.00
            IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
!              For time 18 characters are needed
               FSTR = '(A18)'
               LFIELD = 18
               OUTLIN(LINKAR:LINKAR+LFIELD-1) = CHTIME
            ELSE
               WRITE(NREF,'(A18)') CHTIME
            ENDIF
          ELSE
            IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
              IF (RTYPE.EQ.'TABD') THEN
                FSTR = FLT_TABLE                                          40.13
                LFIELD = FLD_TABLE                                        40.13
              ELSE
                FSTR = '(F13.X)'
                LFIELD = 13
!               NUMDEC is number of decimals in the table for each output quantity
                NUMDEC = MAX (0, 6-NINT(LOG10(ABS(OVHEXP(IVTYPE)))))      40.03
                IF (NUMDEC.GT.9) NUMDEC = 9                               40.00
                WRITE (FSTR(6:6), '(I1)') NUMDEC                          40.00
              ENDIF
              IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 ) THEN
!                    write value into OUTLIN
                WRITE (OUTLIN(LINKAR:LINKAR+LFIELD-1), FMT=FSTR)
     &                 VOQ(IP,VOQR(IVTYPE))
                IF (OVSVTY(IVTYPE).EQ.3) THEN
                  LINKAR = LINKAR + LFIELD + 1
!                       write second component of a vectorial quantity
                  WRITE (OUTLIN(LINKAR:LINKAR+LFIELD-1), FMT=FSTR)
     &                   VOQ(IP,VOQR(IVTYPE)+1)
                ENDIF
              ENDIF
            ELSE
              IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 ) THEN
                 WRITE(NREF,FMT=FLT_SWASH) VOQ(IP,VOQR(IVTYPE))
                 IF (OVSVTY(IVTYPE).EQ.3) THEN
                    WRITE(NREF,FMT=FLT_SWASH) VOQ(IP,VOQR(IVTYPE)+1)
                 ENDIF
              ELSE IF ( IVTYPE.LE.70 ) THEN
                 WRITE(NREF,FMT=FLT_SWASH)
     &                            (VOQK(IP,K,VOQR(IVTYPE)), K = 0, KMAX)
                 IF (OVSVTY(IVTYPE).EQ.3) THEN
                    WRITE(NREF,FMT=FLT_SWASH)
     &                          (VOQK(IP,K,VOQR(IVTYPE)+1), K = 0, KMAX)
                 ENDIF
              ELSE
                 WRITE(NREF,FMT=FLT_SWASH)
     &                            (VOQK(IP,K,VOQR(IVTYPE)), K = 1, KMAX)
                 IF (OVSVTY(IVTYPE).EQ.3) THEN
                    WRITE(NREF,FMT=FLT_SWASH)
     &                          (VOQK(IP,K,VOQR(IVTYPE)+1), K = 1, KMAX)
                 ENDIF
              ENDIF
            ENDIF
          ENDIF
          IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
             IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 ) THEN
                LINKAR = LINKAR + LFIELD + 1
                OUTLIN(LINKAR-1:LINKAR) = '  '
             ENDIF
          ENDIF
  60    CONTINUE
        IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
           IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 )
     &        WRITE (NREF, '(A)') OUTLIN(1:LINKAR-1)
        ENDIF
       END IF                                                             40.51
  70  CONTINUE
!
      RETURN
! * end of subroutine SWTABP *
      END
!MatL4!****************************************************************
!MatL4!
!MatL4      SUBROUTINE SWRMAT ( MROWS , NCOLS, MATNAM, RDATA,
!MatL4     &                    IOUTMA, IREC , IDLA  , DUMVAL )
!MatL4!
!MatL4!****************************************************************
!MatL4!
!MatL4      USE OCPCOMM4                                                        40.41
!MatL4!
!MatL4      IMPLICIT NONE
!MatL4!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
!     Copyright (C) 2010-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!MatL4!
!MatL4!  0. Authors
!MatL4!
!MatL4!     40.30: Marcel Zijlema
!MatL4!     40.41: Marcel Zijlema
!MatL4!
!MatL4!  1. Updates
!MatL4!
!MatL4!     40.30, May 03: New subroutine
!MatL4!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!MatL4!
!MatL4!  2. Purpose
!MatL4!
!MatL4!     Writes block output to a binary file in the MAT-format
!MatL4!     to be used in MATLAB
!MatL4!
!MatL4!  3. Method
!MatL4!
!MatL4!     1) The binary file BINFIL must be opened with the following
!MatL4!        statement:
!MatL4!
!MatL4!        OPEN(UNIT=IOUTMA, FILE=BINFIL, FORM='UNFORMATTED',
!MatL4!             ACCESS='DIRECT', RECL=1)
!MatL4!
!MatL4!        Furthermore, initialize record counter to IREC = 1
!MatL4!
!MatL4!     2) Be sure to close the binary file when there are no more
!MatL4!        matrices to be saved
!MatL4!
!MatL4!     3) The matrix may contain signed infinity and/or Not a Numbers.
!MatL4!        According to the IEEE standard, on a 32-bit machine, the real
!MatL4!        format has an 8-bit biased exponent (=actual exponent increased
!MatL4!        by bias=127) and a 23-bit fraction or mantissa. The leftmost
!MatL4!        bit is the sign bit. Let a fraction, biased exponent and sign
!MatL4!        bit be denoted as F, E and S, respectively. The following
!MatL4!        formats adhere to IEEE standard:
!MatL4!
!MatL4!          S = 0, E = 11111111 and F  = 00 ... 0 : X = +Inf
!MatL4!          S = 1, E = 11111111 and F  = 00 ... 0 : X = -Inf
!MatL4!          S = 0, E = 11111111 and F <> 00 ... 0 : X = NaN
!MatL4!
!MatL4!        Hence, the representation of +Inf equals 2**31 - 2**23. A
!MatL4!        representation of a NaN equals the representation of +Inf
!MatL4!        plus 1.
!MatL4!
!MatL4!        The Cray machine C916 at SARA Information Centre does not
!MatL4!        support the IEEE standard.
!MatL4!
!MatL4!     4) The NaN's or Inf's are indicated by a dummy value as given
!MatL4!        by dumval
!MatL4!
!MatL4!     For more information consult "Appendix - MAT-File Structure"
!MatL4!     of the MATLAB External Data Reference guide (Version 4.2)
!MatL4!
!MatL4!  4. Argument variables
!MatL4!
!MatL4!     DUMVAL      a dummy value meant for indicating NaN
!MatL4!     IDLA        controls lay-out of output (see user manual)
!MatL4!     IOUTMA      unit number of binary MAT-file
!MatL4!     IREC        direct access file record counter
!MatL4!     MATNAM      character array holding the matrix name
!MatL4!     MROWS       a 4-byte integer representing the number of
!MatL4!                 rows in matrix
!MatL4!     NCOLS       a 4-byte integer representing the number of
!MatL4!                 columns in matrix
!MatL4!     RDATA       real array consists of MROWS * NCOLS real
!MatL4!                 elements stored column wise
!MatL4!
!MatL4      INTEGER       MROWS, NCOLS, IDLA, IOUTMA, IREC
!MatL4      REAL          RDATA(*), DUMVAL
!MatL4      CHARACTER*(*) MATNAM
!MatL4!
!MatL4!  5. Parameter variables
!MatL4!
!MatL4!     ---
!MatL4!
!MatL4!  6. Local variables
!MatL4!
!MatL4!     BVAL  :     a byte value
!MatL4!     CHARS :     array to pass character info to MSGERR
!MatL4!     I     :     loop variable
!MatL4!     IENT  :     number of entries
!MatL4!     IF    :     first non-character in string
!MatL4!     IL    :     last non-character in string
!MatL4!     IMAGF :     a 4-byte imaginary flag. Possible values are:
!MatL4!                 0: there is only real data
!MatL4!                 1: the data has also an imaginary part
!MatL4!     IOS   :     auxiliary integer with iostat-number
!MatL4!     ITYPE :     the type flag containing a 4-byte integer whose
!MatL4!                 decimal digits encode storage information.
!MatL4!                 If the integer is represented as ABCD then:
!MatL4!                 "A" indicates the format to write the binary
!MatL4!                 data to a file on the machine. Possible values are:
!MatL4!                   0: Intel based machines (PC 386/486, Pentium)
!MatL4!                   1: Motorola 68000 based machines (Macintosh,
!MatL4!                      HP 9000, SPARC, Apollo, SGI)
!MatL4!                   2: VAX-D format
!MatL4!                   3: VAX-G format
!MatL4!                   4: Cray
!MatL4!                 "B" is always zero
!MatL4!                 "C" indicates which format the data is stored.
!MatL4!                  Possible values are:
!MatL4!                   0: double precision (64 bit) floating point numbers
!MatL4!                   1: single precision (32 bit) floating point numbers
!MatL4!                   2: 32-bit signed integers
!MatL4!                   3: 16-bit signed integers
!MatL4!                   4: 16-bit unsigned integers
!MatL4!                   5: 8-bit unsigned integers
!MatL4!                 "D" indicates the type of data (matrix).
!MatL4!                  Possible values:
!MatL4!                   0: numeric matrix
!MatL4!                   1: textual matrix
!MatL4!                   2: sparse  matrix
!MatL4!     J     :     index
!MatL4!     M     :     loop variable
!MatL4!     MSGSTR:     string to pass message to call MSGERR
!MatL4!     N     :     loop variable
!MatL4!     NAMLEN:     a 4-byte integer representing the number of
!MatL4!                 characters in matrix name plus 1
!MatL4!     NANVAL:     an integer representing Not a Number
!MatL4!
!MatL4      INTEGER I, J, IENT, IF, IL, IOS, M, N
!MatL4      INTEGER BVAL(4), IMAGF, ITYPE, NAMLEN, NANVAL
!MatL4      CHARACTER*20 INTSTR, CHARS
!MatL4      CHARACTER*80 MSGSTR
!MatL4!
!MatL4!  8. Subroutines used
!MatL4!
!MatL4!     INTSTR           Converts integer to string
!MatL4!     MSGERR           Writes error message
!MatL4!     TXPBLA           Removes leading and trailing blanks in string
!MatL4!     SWI2B            Calculates 32-bit representation of an
!MatL4!                      integer number
!MatL4!     SWR2B            Calculates 32-bit representation of a
!MatL4!                      floating-point number
!MatL4!
!MatL4!  9. Subroutines calling
!MatL4!
!MatL4!     ---
!MatL4!
!MatL4! 10. Error messages
!MatL4!
!MatL4!     ---
!MatL4!
!MatL4! 11. Remarks
!MatL4!
!MatL4!     ---
!MatL4!
!MatL4! 12. Structure
!MatL4!
!MatL4!     set Not a Number
!MatL4!
!MatL4!     set some flags
!MatL4!
!MatL4!     write header consisting of ITYPE, MROWS, NCOLS, IMAGF, NAMLEN and
!MatL4!     name of matrix MATNAM
!MatL4!
!MatL4!     write matrix
!MatL4!
!MatL4!     if necessary, give message that error occurred while writing file
!MatL4!
!MatL4! 13. Source text
!MatL4!
!MatL4      SAVE IENT
!MatL4      DATA IENT/0/
!MatL4      IF (LTRACE) CALL STRACE (IENT,'SWRMAT')
!MatL4
!MatL4!     --- set Not a Number
!MatL4
!MatL4      NANVAL = 255 * 2**23 + 1
!MatL4
!MatL4!     --- set some flags
!MatL4
!MatL4      ITYPE = 1010
!MatL4      IMAGF = 0
!MatL4      IOS   = 0
!MatL4
!MatL4!     --- write header consisting of ITYPE, MROWS, NCOLS, IMAGF,
!MatL4!         NAMLEN and name of matrix MATNAM
!MatL4!         the name should be ended by zero-byte terminator
!MatL4
!MatL4      CALL SWI2B ( ITYPE, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      CALL SWI2B ( MROWS, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      CALL SWI2B ( NCOLS, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      CALL SWI2B ( IMAGF, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      CALL TXPBLA(MATNAM,IF,IL)
!MatL4      NAMLEN = IL - IF + 2
!MatL4      CALL SWI2B ( NAMLEN, BVAL )
!MatL4      DO I = 1, 4
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4
!MatL4      DO I = IF, IL
!MatL4         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) MATNAM(I:I)
!MatL4         IREC = IREC + 1
!MatL4      END DO
!MatL4      IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(0)
!MatL4      IREC = IREC + 1
!MatL4
!MatL4!     --- write matrix
!MatL4
!MatL4      DO M = 1, NCOLS
!MatL4         DO N = 1, MROWS
!MatL4            IF ( IDLA.EQ.1 ) THEN
!MatL4               J = (MROWS-N)*NCOLS + M
!MatL4            ELSE
!MatL4               J = (N-1)*NCOLS + M
!MatL4            END IF
!MatL4            IF (RDATA(J).NE.DUMVAL ) THEN
!MatL4               CALL SWR2B ( RDATA(J), BVAL )
!MatL4            ELSE IF (.NOT. DUMVAL.NE.0. ) THEN
!MatL4               CALL SWR2B ( RDATA(J), BVAL )
!MatL4            ELSE
!MatL4               CALL SWI2B ( NANVAL, BVAL )
!MatL4            END IF
!MatL4            DO I = 1, 4
!MatL4              IF (IOS.EQ.0)
!MatL4     &                   WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
!MatL4              IREC = IREC + 1
!MatL4            END DO
!MatL4         END DO
!MatL4      END DO
!MatL4
!MatL4!     --- if necessary, give message that error occurred while writing file
!MatL4
!MatL4      IF ( IOS.NE.0 ) THEN
!MatL4         CHARS = INTSTR(IOS)
!MatL4         CALL TXPBLA(CHARS,IF,IL)
!MatL4         MSGSTR = 'Error while writing binary MAT-file - '//
!MatL4     &            'IOSTAT number is '//CHARS(IF:IL)
!MatL4         CALL MSGERR ( 4, MSGSTR )
!MatL4         RETURN
!MatL4      END IF
!MatL4
!MatL4      RETURN
!MatL4      END
!****************************************************************
!
      SUBROUTINE SWRMAT ( MROWS , NCOLS, MATNAM, RDATA,
     &                    IOUTMA, IREC , IDLA  , DUMVAL )
!
!****************************************************************
!
      USE OCPCOMM2                                                        41.08
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
!     Copyright (C) 2010-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     41.08: Pieter Smit
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.08, Aug. 09: adapted to write binary file in the Level 5 MAT-file format
!
!  2. Purpose
!
!     Writes block output to a binary file in the MAT-format
!     to be used in MATLAB
!
!  3. Method
!
!     1) The binary file BINFIL must be opened with the following
!        statement:
!
!        OPEN(UNIT=IOUTMA, FILE=BINFIL, FORM='UNFORMATTED',
!             ACCESS='DIRECT', RECL=4)
!
!        Furthermore, initialize record counter to IREC = 1
!
!     2) Be sure to close the binary file when there are no more
!        matrices to be saved
!
!     3) The matrix may contain signed infinity and/or Not a Numbers.
!        According to the IEEE 754 standard, on a 32-bit machine, the real
!        format has an 8-bit biased exponent (=actual exponent increased
!        by bias=127) and a 23-bit fraction or mantissa. The leftmost
!        bit is the sign bit. Let a fraction, biased exponent and sign
!        bit be denoted as F, E and S, respectively. The following
!        formats adhere to IEEE standard:
!
!          S = 0, E = 11111111 and F  = 00 ... 0 : X = +Inf
!          S = 1, E = 11111111 and F  = 00 ... 0 : X = -Inf
!          S = 0, E = 11111111 and F <> 00 ... 0 : X = NaN
!
!        Hence, the representation of +Inf equals 2**31 - 2**23. A
!        representation of a NaN equals the representation of +Inf
!        plus 1.
!
!     4) The NaN's or Inf's are indicated by a dummy value as given
!        by dumval
!
!     For more information on the Level 5 MAT-file format consult
!     document "MAT-File Format" of MathWorks
!
!  4. Argument variables
!
!     DUMVAL      a dummy value meant for indicating NaN
!     IDLA        controls lay-out of output (see user manual)
!     IOUTMA      unit number of binary MAT-file
!     IREC        direct access file record counter
!     MATNAM      character array holding the matrix name
!     MROWS       a 4-byte integer representing the number of
!                 rows in matrix
!     NCOLS       a 4-byte integer representing the number of
!                 columns in matrix
!     RDATA       real array consists of MROWS * NCOLS real
!                 elements stored column wise
!
      INTEGER       MROWS, NCOLS, IDLA, IOUTMA, IREC
      REAL          RDATA(*), DUMVAL
      CHARACTER*(*) MATNAM
!
!  5. Parameter variables
!
!     BlockSize   size of matlab data segment
!     DataSize    number of bytes written per write statement
!     HeaderSize  size of the header in bytes
!     mChar       character data
!     mInt32      signed   integer*4
!     mUInt32     unsigned integer*8
!     mSingle     real
!
!     --- standard sizes
!
      INTEGER, PARAMETER :: DataSize   = 4
      INTEGER, PARAMETER :: HeaderSize = 128
      INTEGER, PARAMETER :: BlockSize  = 8
!
!     --- Matlab data types
!
      INTEGER, PARAMETER :: mChar      = 1
      INTEGER, PARAMETER :: mInt32     = 5
      INTEGER, PARAMETER :: mUInt32    = 6
      INTEGER, PARAMETER :: mSingle    = 7
!
!  6. Local variables
!
!     CTMP  :     a temporary character array
!     HEADER:     header of binary MAT-file
!     I     :     loop variable
!     IENT  :     number of entries
!     IOS   :     auxiliary integer with iostat-number
!     IRECS :     size of array including tags and flags
!     J     :     index
!     M     :     loop variable
!     MSGSTR:     string to pass message to call MSGERR
!     N     :     loop variable
!     NAMLEN:     a 4-byte integer representing the number of
!                 characters in matrix name
!     NANVAL:     an integer representing Not a Number
!     NTOT  :     size of data array
!
      INTEGER I, J, IENT, IOS, M, N, NTOT
      INTEGER NAMLEN, NANVAL
      INTEGER, SAVE :: IRECS
      CHARACTER*80 MSGSTR
      CHARACTER(LEN=HeaderSize) HEADER
      CHARACTER(LEN=BlockSize) CTMP

!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     set Not a Number
!     length of name matrix
!     size of data array
!     write header once
!     array name
!     write matrix
!     write the size of the array
!     if necessary, give message that error occurred while writing file
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRMAT')

      IOS = 0

!     --- set Not a Number

      NANVAL = 255 * 2**23 + 1

!     --- length of name matrix

      NAMLEN = LEN_TRIM(MATNAM)

!     --- size of data array

      NTOT = MROWS * NCOLS

!     --- descriptive header

      WRITE (HEADER, '(6A)') 'Data produced by SWASH version ',
     &                       TRIM(VERTXT),'; project: ',TRIM(PROJID),
     &                       '; run number: ',PROJNR

!     --- data offset

      HEADER(117:124) = CHAR(ICHAR(' '))

!     --- version

      HEADER(125:126) = CHAR(0) // CHAR(1)

!     --- endian indicator

      WRITE(HEADER(127:128),'(A)') INT(19785,KIND=2)

!     --- write header once

      IF ( IREC.EQ.1 ) THEN
         DO I = 1, HeaderSize/DataSize
            J = DataSize*(I-1) + 1
            IF ( IOS.EQ.0 )
     &         WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) HEADER(J:J+DataSize-1)
            IREC = IREC + 1
         END DO
      END IF

!     --- array tag

      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 14                ! Matlab array
      IREC = IREC + 1
      IRECS = IREC
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
      IREC = IREC + 1

!     --- array flags

      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mInt32
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 2*DataSize
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 7                 ! single precision array
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
      IREC = IREC + 1

      IF ( MOD(2,BlockSize/DataSize).NE.0 ) THEN
         IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
         IREC = IREC + 1
      END IF

!     --- dimensions array

      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mInt32
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 2*DataSize
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) MROWS
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) NCOLS
      IREC = IREC + 1

      IF ( MOD(2,BlockSize/DataSize).NE.0 ) THEN
         IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
         IREC = IREC + 1
      END IF

!     --- array name

      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mChar
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) NAMLEN
      IREC = IREC + 1

      I = 1
      DO
        CTMP(1:8) = CHAR(ICHAR(' '))
        IF ( I.GT.NAMLEN ) THEN
           EXIT
        ELSE IF ( I+BlockSize.LE.NAMLEN ) THEN
           CTMP(1:8) = MATNAM(I:I+BlockSize-1)
        ELSE
           CTMP(1:NAMLEN-I+1) = MATNAM(I:NAMLEN)
        END IF

        IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CTMP(1:4)
        IREC = IREC + 1
        IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CTMP(5:8)
        I = I + BlockSize
        IREC = IREC + 1
      END DO

!     --- write matrix

      IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) mSingle
      IREC = IREC + 1
      IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) NTOT*DataSize
      IREC = IREC + 1

      DO M = 1, NCOLS
         DO N = 1, MROWS
            IF ( IDLA.EQ.1 ) THEN
               J = (MROWS-N)*NCOLS + M
            ELSE
               J = (N-1)*NCOLS + M
            END IF
            IF (RDATA(J).NE.DUMVAL ) THEN
               WRITE (IOUTMA,REC=IREC) RDATA(J)
            ELSE IF (.NOT. DUMVAL.NE.0. ) THEN
               WRITE (IOUTMA,REC=IREC) RDATA(J)
            ELSE
               WRITE (IOUTMA,REC=IREC) NANVAL
            END IF
            IREC = IREC + 1
         END DO
      END DO

      IF ( MOD(NTOT,BlockSize/DataSize).NE.0 ) THEN
         IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) 0.
         IREC = IREC + 1
      END IF

!     --- write the size of the array

      IF ( IOS.EQ.0 )
     &   WRITE (IOUTMA,REC=IRECS,IOSTAT=IOS) (IREC-IRECS-1)*DataSize

!     --- if necessary, give message that error occurred while writing file

      IF ( IOS.NE.0 ) THEN
         WRITE (MSGSTR, '(A,I5)')
     &            'Error while writing binary MAT-file - '//
     &            'IOSTAT number is ', IOS
         CALL MSGERR( 4, TRIM(MSGSTR) )
         RETURN
      END IF

      RETURN
      END
