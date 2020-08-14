!
!     SWAN - routines for distributed-memory approach based on MPI
!
!  Contents of this file
!
!     SWINITMPI
!     SWEXITMPI
!     SWSYNC
!     SWSENDNB
!     SWSENDNBD
!     SWRECVNB
!     SWRECVNBD
!     SWBROADC
!     SWBROADCD
!     SWGATHER
!     SWGATHERD
!     SWREDUCE
!     SWREDUCI
!     SWREDUCR
!     SWREDUCD
!     SWSTRIP
!     SWORB
!     SWPARTIT
!     SWBLADM
!     SWCORN
!     SWDECOMP
!     SWEXCHG
!     SWEXCHGI
!     SWEXCHGD
!     SWRECVFLD
!     SWSENDFLD
!     SWCOLLECT
!     SWCOLOUT
!     SWCOLTAB
!     SWCOLBLK
!
!****************************************************************
!
      SUBROUTINE SWINITMPI
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Join parallel application
!
!  3. Method
!
!     Start MPI and initialize some variables
!
!  4. Argument variables
!
!     ---
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_COMM_RANK    Get rank of processes in MPI communication context
!MPI!     MPI_COMM_SIZE    Get number of processes in MPI communication context
!MPI!     MPI_INIT         Enroll in MPI
!     MSGERR           Writes error message
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     Main program SWAN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Start MPI and initialize some common variables in module M_PARALL
!
! 13. Source text
!
      LEVERR = 0
      MAXERR = 1
      ITRACE = 0

!MPI!     --- enroll in MPI
!MPI
!MPI      CALL MPI_INIT ( IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- initialize common variables

      INODE = 0
      NPROC = 1

!MPI!     --- get node number INODE
!MPI
!MPI      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, INODE, IERR )
      INODE = INODE + 1
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
!MPI
!MPI!     --- determine total number of processes
!MPI
!MPI      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, NPROC, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- determine whether this is a parallel run or not

      IF ( NPROC.GT.1 ) THEN
         PARLL = .TRUE.
      ELSE
         PARLL = .FALSE.
      END IF

!MPI!     --- define MPI constants for communication within SWAN
!MPI
!MPI      SWINT  = MPI_INTEGER
!MPI      SWREAL = MPI_REAL
!MPI      SWDBLE = MPI_DOUBLE_PRECISION
!MPI      SWMAX  = MPI_MAX
!MPI      SWMIN  = MPI_MIN
!MPI      SWSUM  = MPI_SUM

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXITMPI
!
!****************************************************************
!
!MPI      USE MPI
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
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Exit parallel application
!
!  3. Method
!
!     Wrapper for MPI_FINALIZE
!
!  4. Argument variables
!
!     ---
!
!  6. Local variables
!
!     IERR    :   error value of MPI call
!     PARALMPI:   if true, parallel process is carried out with MPI
!
      INTEGER IERR
      LOGICAL PARALMPI
!
!  8. Subroutines used
!
!MPI!     MPI_ABORT        Abort MPI if severe error occurs
!MPI!     MPI_BARRIER      Blocks until all nodes have called this routine
!MPI!     MPI_INITIALIZED  Indicates whether MPI_Init has been called
!MPI!     MPI_FINALIZE     Cleans up the MPI state and exits
!
!  9. Subroutines calling
!
!     Main program SWAN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     if MPI has been initialized
!MPI!        synchronize nodes
!MPI!        if severe error
!MPI!           abort MPI
!MPI!        else
!MPI!           close MPI
!MPI!
! 13. Source text
!
!MPI      CALL MPI_INITIALIZED ( PARALMPI, IERR )
!MPI      IF ( PARALMPI ) THEN
!MPI
!MPI         CALL MPI_BARRIER ( MPI_COMM_WORLD, IERR )
!MPI
!MPI         IF ( LEVERR.GE.4 ) THEN
!MPI
!MPI!        --- in case of a severe error abort all MPI processes
!MPI
!MPI            CALL MPI_ABORT ( MPI_COMM_WORLD, LEVERR, IERR )
!MPI
!MPI         ELSE
!MPI
!MPI!        --- otherwise stop MPI operations on this computer
!MPI
!MPI            CALL MPI_FINALIZE ( IERR )
!MPI
!MPI         END IF
!MPI
!MPI      END IF
!MPI
      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSYNC
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Jan. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Synchronize nodes
!
!  3. Method
!
!     Wrapper for MPI_BARRIER
!
!  4. Argument variables
!
!     ---
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_BARRIER      Blocks until all nodes have called this routine
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     Blocks until all nodes have called MPI_BARRIER routine.
!     In this way, all nodes are synchronized
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSYNC')

!MPI!     --- blocks until all nodes have called this routine
!MPI
!MPI      CALL MPI_BARRIER ( MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDNB ( IPTR, ILEN, ITYPE, IDEST, ITAG )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Data is sent to a neighbour
!
!  3. Method
!
!     Wrapper for MPI_SEND
!
!  4. Argument variables
!
!     IDEST       rank of the destination process
!     ILEN        length of array to be sent
!     IPTR        pointer to first element of array to be sent
!     ITAG        message type
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE, IDEST, ITAG
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_SEND         Immediately sends the data in the active
!MPI!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHG
!     SWEXCHGI
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     Data is sent to a neighbour with command MPI_SEND
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDNB')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!MPI      CALL MPI_SEND ( IPTR, ILEN, ITYPE, IDEST-1,
!MPI     &                ITAG, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDNBD ( ARR, ILEN, IDEST, ITAG )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
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
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 17: New subroutine
!
!  2. Purpose
!
!     Data of double precision is sent to a neighbour
!
!  3. Method
!
!     Wrapper for MPI_SEND
!
!  4. Argument variables
!
!     ARR         double precision array
!     IDEST       rank of the destination process
!     ILEN        length of array to be sent
!     ITAG        message type
!
      INTEGER ILEN, IDEST, ITAG
      REAL*8  ARR(ILEN)
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_SEND         Immediately sends the data in the active
!MPI!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHGD
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     Data is sent to a neighbour with command MPI_SEND
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDNBD')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!MPI      CALL MPI_SEND ( ARR, ILEN, SWDBLE, IDEST-1,
!MPI     &                ITAG, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVNB ( IPTR, ILEN, ITYPE, ISOURCE, ITAG )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Data is received from a neighbour
!
!  3. Method
!
!     Wrapper for MPI_RECV
!
!  4. Argument variables
!
!     ILEN        length of array to be received
!     IPTR        pointer to first element of array to be received
!     ISOURCE     rank of the source process
!     ITAG        message type
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE, ISOURCE, ITAG
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!MPI!     ISTAT :     MPI status array
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
!MPI      INTEGER      ISTAT(MPI_STATUS_SIZE)
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_RECV         Immediately receives the data in the active
!MPI!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHG
!     SWEXCHGI
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
!MPI!     Data is received from a neighbour with command MPI_RECV
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVNB')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!MPI      CALL MPI_RECV ( IPTR, ILEN, ITYPE, ISOURCE-1, ITAG,
!MPI     &                MPI_COMM_WORLD, ISTAT, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVNBD ( ARR, ILEN, ISOURCE, ITAG )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
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
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 17: New subroutine
!
!  2. Purpose
!
!     Data of double precision is received from a neighbour
!
!  3. Method
!
!     Wrapper for MPI_RECV
!
!  4. Argument variables
!
!     ARR         double precision array
!     ILEN        length of array to be received
!     ISOURCE     rank of the source process
!     ITAG        message type
!
      INTEGER ILEN, ISOURCE, ITAG
      REAL*8  ARR(ILEN)
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF1   :     first non-character in string1
!     IF2   :     first non-character in string2
!     IL1   :     last non-character in string1
!     IL2   :     last non-character in string2
!MPI!     ISTAT :     MPI status array
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF1, IF2, IL1, IL2
!MPI      INTEGER      ISTAT(MPI_STATUS_SIZE)
      CHARACTER*20 INTSTR, CHARS(2)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_RECV         Immediately receives the data in the active
!MPI!                      MPI message buffer
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWEXCHGD
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
!MPI!     Data is received from a neighbour with command MPI_RECV
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVNBD')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!MPI      CALL MPI_RECV ( ARR, ILEN, SWDBLE, ISOURCE-1, ITAG,
!MPI     &                MPI_COMM_WORLD, ISTAT, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS(1) = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS(1),IF1,IL1)
!MPI         CHARS(2) = INTSTR(INODE)
!MPI         CALL TXPBLA(CHARS(2),IF2,IL2)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(1)(IF1:IL1)//
!MPI     &            ' and node number is '//CHARS(2)(IF2:IL2)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBROADC ( IPTR, ILEN, ITYPE )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Broadcasts data from the master to all other processes
!
!  3. Method
!
!     Wrapper for MPI_BCAST
!
!  4. Argument variables
!
!     ILEN        length of array to be sent
!     IPTR        pointer to first element of array to be sent
!     ITYPE       type of data
!
      INTEGER IPTR, ILEN, ITYPE
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_BCAST        Broadcasts a message from the master
!MPI!                      to all other processes of the group
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SNEXTI
!     FLFILE
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     Broadcasts data from the master to all other nodes
!MPI!     with command MPI_BCAST
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBROADC')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(201)
!MPI      CALL MPI_BCAST ( IPTR, ILEN, ITYPE, MASTER-1,
!MPI     &                 MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
!TIMG      CALL SWTSTO(201)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBROADCD ( ARR, ILEN )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
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
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 17: New subroutine
!
!  2. Purpose
!
!     Broadcasts data of double precision from the master to all other processes
!
!  3. Method
!
!     Wrapper for MPI_BCAST
!
!  4. Argument variables
!
!     ARR         double precision array
!     ILEN        length of array to be sent
!
      INTEGER ILEN
      REAL*8  ARR(ILEN)
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_BCAST        Broadcasts a message from the master
!MPI!                      to all other processes of the group
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     DAC2
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!MPI!     Broadcasts double precisions from the master to all other nodes
!MPI!     with command MPI_BCAST
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBROADCD')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(201)
!MPI      CALL MPI_BCAST ( ARR, ILEN, SWDBLE, MASTER-1,
!MPI     &                 MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
!TIMG      CALL SWTSTO(201)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWGATHER ( IOPTR, IOLEN, IIPTR, IILEN, ITYPE )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Gathers different amounts of data from each processor
!     to the master
!
!  3. Method
!
!     Wrapper for MPI_GATHERV
!
!  4. Argument variables
!
!     IILEN       length of input array
!     IIPTR       pointer to first element of input array (local)
!     IOLEN       length of output array
!     IOPTR       pointer to first element of output array (global)
!     ITYPE       type of data
!
      INTEGER IILEN, IIPTR, IOLEN, IOPTR, ITYPE
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     ICOUNT:     array specifying array size of data received
!                 from each processor
!     IDSPLC:     array specifying the starting address of the
!                 incoming data from each processor, relative
!                 to the global array
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER                 I, IENT, IERR, IF, IL
      INTEGER, ALLOCATABLE :: ICOUNT(:), IDSPLC(:)
      CHARACTER*20            INTSTR, CHARS
      CHARACTER*80            MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_GATHER       Gathers data from all nodes to the master
!MPI!     MPI_GATHERV      Gathers different amounts of data from
!MPI!                      all nodes to the master
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!MPI!     gather the array sizes to the master
!MPI!
!     check whether enough space has been allocated
!     for gathered data
!
!     calculate starting address of each local array
!     with respect to the global array
!
!MPI!     gather different amounts of data from each processor
!MPI!     to the master
!MPI!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWGATHER')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IF (INODE.EQ.MASTER) THEN
         ALLOCATE(ICOUNT(0:NPROC-1))
         ALLOCATE(IDSPLC(0:NPROC-1))
      END IF

!MPI!     --- gather the array sizes to the master
!MPI
!MPI      CALL MPI_GATHER( IILEN, 1, SWINT, ICOUNT, 1, SWINT,
!MPI     &                 MASTER-1, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- check whether enough space has been allocated
!         for gathered data

      IF (INODE.EQ.MASTER) THEN
         IF ( SUM(ICOUNT).GT.IOLEN ) THEN
            CALL MSGERR(4,
     &                  'Not enough space allocated for gathered data')
            RETURN
         END IF
      END IF

!     --- calculate starting address of each local array
!         with respect to the global array

      IF (INODE.EQ.MASTER) THEN
         IDSPLC(0) = 0
         DO I = 1, NPROC-1
            IDSPLC(I) = ICOUNT(I-1) + IDSPLC(I-1)
         END DO
      END IF

!MPI!     --- gather different amounts of data from each processor
!MPI!         to the master
!MPI
!MPI      CALL MPI_GATHERV( IIPTR, IILEN, ITYPE, IOPTR, ICOUNT, IDSPLC,
!MPI     &                  ITYPE, MASTER-1, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      IF (INODE.EQ.MASTER) DEALLOCATE(ICOUNT,IDSPLC)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWGATHERD ( ARRO, IOLEN, ARRI, IILEN )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
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
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 17: New subroutine
!
!  2. Purpose
!
!     Gathers different amounts of data of double precision from each processor
!     to the master
!
!  3. Method
!
!     Wrapper for MPI_GATHERV
!
!  4. Argument variables
!
!     ARRI        double precision input array (local)
!     ARRO        double precision output array (global)
!     IILEN       length of input array
!     IOLEN       length of output array
!
      INTEGER IILEN, IOLEN
      REAL*8  ARRI(IILEN), ARRO(IOLEN)
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     ICOUNT:     array specifying array size of data received
!                 from each processor
!     IDSPLC:     array specifying the starting address of the
!                 incoming data from each processor, relative
!                 to the global array
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER                 I, IENT, IERR, IF, IL
      INTEGER, ALLOCATABLE :: ICOUNT(:), IDSPLC(:)
      CHARACTER*20            INTSTR, CHARS
      CHARACTER*80            MSGSTR
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_GATHER       Gathers data from all nodes to the master
!MPI!     MPI_GATHERV      Gathers different amounts of data from
!MPI!                      all nodes to the master
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     DAC2
!
! 10. Error messages
!
!     ---
!
! 12. Structure
!
!     if not parallel, return
!
!MPI!     gather the array sizes to the master
!MPI!
!     check whether enough space has been allocated
!     for gathered data
!
!     calculate starting address of each local array
!     with respect to the global array
!
!MPI!     gather different amounts of data from each processor
!MPI!     to the master
!MPI!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWGATHERD')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IF (INODE.EQ.MASTER) THEN
         ALLOCATE(ICOUNT(0:NPROC-1))
         ALLOCATE(IDSPLC(0:NPROC-1))
      END IF

!MPI!     --- gather the array sizes to the master
!MPI
!MPI      CALL MPI_GATHER( IILEN, 1, SWINT, ICOUNT, 1, SWINT,
!MPI     &                 MASTER-1, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

!     --- check whether enough space has been allocated
!         for gathered data

      IF (INODE.EQ.MASTER) THEN
         IF ( SUM(ICOUNT).GT.IOLEN ) THEN
            CALL MSGERR(4,
     &                  'Not enough space allocated for gathered data')
            RETURN
         END IF
      END IF

!     --- calculate starting address of each local array
!         with respect to the global array

      IF (INODE.EQ.MASTER) THEN
         IDSPLC(0) = 0
         DO I = 1, NPROC-1
            IDSPLC(I) = ICOUNT(I-1) + IDSPLC(I-1)
         END DO
      END IF

!MPI!     --- gather different amounts of data from each processor
!MPI!         to the master
!MPI
!MPI      CALL MPI_GATHERV( ARRI, IILEN, SWDBLE, ARRO, ICOUNT, IDSPLC,
!MPI     &                  SWDBLE, MASTER-1, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF

      IF (INODE.EQ.MASTER) DEALLOCATE(ICOUNT,IDSPLC)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCE ( IPTR, ILEN, ITYPE, ITYPRD )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.96, Dec. 08: call SWREDUCI/R instead of passing startaddress of the array
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on
!     array (IPTR) of type ITYPE to collect values from
!     all processes
!
!  4. Argument variables
!
!     ILEN        length of array to be collect
!     IPTR        pointer to first element of array to be collect
!     ITYPE       type of data
!     ITYPRD      type of reduction
!
      INTEGER IPTR, ILEN, ITYPE, ITYPRD
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     IENT  :     number of entries
!
      INTEGER IENT
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWCOMP
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
!     Performs a global reduction of data across all nodes
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCE')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN
!
!     --- actual reduction of field array based on its type
      IF ( ITYPE.EQ.SWINT ) THEN
         CALL SWREDUCI ( IPTR, ILEN, ITYPRD )
      ELSE IF ( ITYPE.EQ.SWREAL ) THEN
         CALL SWREDUCR ( IPTR, ILEN, ITYPRD )
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCI ( IARR, ILEN, ITYPRD )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
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
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 08: New subroutine
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on integer
!     array IARR to collect values from all processes
!
!  3. Method
!
!     Wrapper for MPI_ALLREDUCE
!
!  4. Argument variables
!
!     IARR        integer array
!     ILEN        length of array to be collect
!     ITYPRD      type of reduction
!
      INTEGER ILEN, ITYPRD
      INTEGER IARR(ILEN)
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ITEMP :     temporary array to store collected data
!     MSGSTR:     string to pass message to call MSGERR
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR

      INTEGER, ALLOCATABLE :: ITEMP(:)
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_ALLREDUCE    Combines values from all processes and
!MPI!                      distribute the result back to all processes
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWREDUCE
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
!MPI!     Performs a global reduction of integers across all nodes
!MPI!     with command MPI_ALLREDUCE
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCI')

      ALLOCATE(ITEMP(ILEN))

!TIMG      CALL SWTSTA(202)
!MPI      CALL MPI_ALLREDUCE ( IARR, ITEMP, ILEN, SWINT,
!MPI     &                     ITYPRD, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
      IARR = ITEMP
!TIMG      CALL SWTSTO(202)

      DEALLOCATE(ITEMP)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCR ( ARR, ILEN, ITYPRD )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
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
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 08: New subroutine
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on real
!     array ARR to collect values from all processes
!
!  3. Method
!
!     Wrapper for MPI_ALLREDUCE
!
!  4. Argument variables
!
!     ARR         real array
!     ILEN        length of array to be collect
!     ITYPRD      type of reduction
!
      INTEGER ILEN, ITYPRD
      REAL    ARR(ILEN)
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!     TEMP  :     temporary array to store collected data
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR

      REAL, ALLOCATABLE :: TEMP(:)
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_ALLREDUCE    Combines values from all processes and
!MPI!                      distribute the result back to all processes
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWREDUCE
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
!MPI!     Performs a global reduction of reals across all nodes
!MPI!     with command MPI_ALLREDUCE
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCR')

      ALLOCATE(TEMP(ILEN))

!TIMG      CALL SWTSTA(202)
!MPI      CALL MPI_ALLREDUCE ( ARR, TEMP, ILEN, SWREAL,
!MPI     &                     ITYPRD, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
      ARR = TEMP
!TIMG      CALL SWTSTO(202)

      DEALLOCATE(TEMP)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWREDUCD ( ARR, ILEN, ITYPRD )
!
!****************************************************************
!
!MPI      USE MPI
      USE OCPCOMM4
      USE M_PARALL
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
!     40.96: Marcel Zijlema
!
!  1. Updates
!
!     40.96, Dec. 08: New subroutine
!
!  2. Purpose
!
!     Performs a global reduction of type ITYPRD on double precision
!     array ARR to collect values from all processes
!
!  3. Method
!
!     Wrapper for MPI_ALLREDUCE
!
!  4. Argument variables
!
!     ARR         double precision array
!     ILEN        length of array to be collect
!     ITYPRD      type of reduction
!
      INTEGER ILEN, ITYPRD
      REAL*8  ARR(ILEN)
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     IENT  :     number of entries
!     IERR  :     error value of MPI call
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     MSGSTR:     string to pass message to call MSGERR
!     TEMP  :     temporary array to store collected data
!
      INTEGER      IENT, IERR, IF, IL
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR

      REAL*8, ALLOCATABLE :: TEMP(:)
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!MPI!     MPI_ALLREDUCE    Combines values from all processes and
!MPI!                      distribute the result back to all processes
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWREDUCE
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
!MPI!     Performs a global reduction of double precisions across all nodes
!MPI!     with command MPI_ALLREDUCE
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWREDUCD')

      ALLOCATE(TEMP(ILEN))

!TIMG      CALL SWTSTA(202)
!MPI      CALL MPI_ALLREDUCE ( ARR, TEMP, ILEN, SWDBLE,
!MPI     &                     ITYPRD, MPI_COMM_WORLD, IERR )
!MPI      IF ( IERR.NE.MPI_SUCCESS ) THEN
!MPI         CHARS = INTSTR(IERR)
!MPI         CALL TXPBLA(CHARS,IF,IL)
!MPI         MSGSTR = 'MPI produces some internal error - '//
!MPI     &            'return code is '//CHARS(IF:IL)
!MPI         CALL MSGERR ( 4, MSGSTR )
!MPI         RETURN
!MPI      END IF
      ARR = TEMP
!TIMG      CALL SWTSTO(202)

      DEALLOCATE(TEMP)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSTRIP ( IPOWN, IDIR, IPART, NPART, LPARTS,
     &                     MXC  , MYC )
!
!****************************************************************
!
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Performs a stripwise partitioning with straight interfaces
!
!  3. Method
!
!     This method is described in Ph.D. Thesis of M. Roest
!     entitled:
!     Partitioning for parallel finite difference computations
!     in coastal water simulation, DUT, 1997
!
!  4. Argument variables
!
!     IDIR        direction of cutting
!                 1 = row
!                 2 = column
!     IPART       part number that must be partitioned
!     IPOWN       array giving the subdomain number of each gridpoint
!     LPARTS      list of parts to be created
!                    lparts(1,i) = number of i-th part to be created
!                    lparts(2,i) = size of i-th part to be created
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!     NPART       number of parts to be created
!
      INTEGER   IDIR, IPART, MXC, MYC, NPART
      INTEGER   IPOWN(*)
      INTEGER*8 LPARTS(2,*)
!
!  6. Local variables
!
!     IC    :     index of (IX,IY)-point
!     ICC   :     index of (IX,IY)-point
!     IENT  :     number of entries
!     INCX  :     increment for adressing: 1 for x-dir, MXC for y-dir
!     INCY  :     increment for adressing: MXC for x-dir, 1 for y-dir
!     IX    :     index in x-direction
!     IY    :     index in y-direction
!     IYY   :     index in y-direction
!     JPART :     a part counter
!     MXCI  :     maximum counter of gridpoints in x/y-direction
!     MYCI  :     maximum counter of gridpoints in y/x-direction
!     NBACK :     number of points in a row already assigned to a part
!     NFORW :     number of points in a row remaining to be assigned
!     NINPRT:     number of points currently assigned to a new part
!
      INTEGER IC, ICC, IENT, INCX, INCY, IX, IY, IYY, JPART,
     &        MXCI, MYCI, NBACK, NFORW, NINPRT
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWORB
!     SWPARTIT
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
!     depending on IDIR, determine indirect addressing in IPOWN
!     start by creating the first part, which is currently empty
!     for all points in IPOWN do
!        if point belongs to the part that must be partitioned then
!           when current part has reached its planned size
!              see how many points in this row remain to be assigned
!              if no more points to be assigned, go on to next part
!              else, if majority of row has been assigned take the rest
!              else, leave this row to next part
!           assign point IC to part that is currently being created
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSTRIP')

!     --- depending on IDIR, determine indirect addressing in IPOWN

      IF ( IDIR.EQ.1 ) THEN
         MXCI = MYC
         MYCI = MXC
         INCX = MXC
         INCY = 1
      ELSE IF ( IDIR.EQ.2 ) THEN
         MXCI = MXC
         MYCI = MYC
         INCX = 1
         INCY = MXC
      END IF

!     --- start by creating the first part, which is currently empty

      JPART  = 1
      NINPRT = 0

!     --- for all points in IPOWN do

      DO IX = 1, MXCI
         DO IY = 1, MYCI

            IC = IX*INCX + IY*INCY - MXC

!           --- if this point belongs to the part that must be partitioned then

            IF ( IPOWN(IC).EQ.IPART ) THEN

!              --- when current part has reached its planned size

               IF ( NINPRT.GE.LPARTS(2,JPART) ) THEN

!                 --- see how many points in this row have been assigned

                  NBACK = 0
                  DO IYY = 1, IY-1

                     ICC = IX*INCX + IYY*INCY - MXC
                     IF (IPOWN(ICC).EQ.LPARTS(1,JPART)) NBACK = NBACK +1

                  END DO

!                 --- see how many points in this row remain to be assigned

                  NFORW = 0
                  DO IYY = IY, MYCI

                     ICC = IX*INCX + IYY*INCY - MXC
                     IF (IPOWN(ICC).EQ.IPART) NFORW = NFORW +1

                  END DO

!                 --- if no more points to be assigned, go on to next part

                  IF ( NFORW.EQ.0 ) THEN

                     JPART  = JPART + 1
                     NINPRT = 0

                  ELSE IF ( (NBACK-NFORW).GT.0 ) THEN

!                    --- if majority of row has been assigned take the rest

                     LPARTS(2,JPART  ) = LPARTS(2,JPART  ) + NFORW
                     LPARTS(2,JPART+1) = LPARTS(2,JPART+1) - NFORW

                  ELSE
!                    --- else, leave this row to next part

                     LPARTS(2,JPART  ) = LPARTS(2,JPART  ) - NBACK
                     LPARTS(2,JPART+1) = LPARTS(2,JPART+1) + NBACK

                     DO IYY = 1, IY-1

                        ICC = IX*INCX + IYY*INCY - MXC
                        IF ( IPOWN(ICC).EQ.LPARTS(1,JPART) ) THEN
                           IPOWN(ICC) = LPARTS(1,JPART+1)
                        END IF

                     END DO

                     JPART  = JPART + 1
                     NINPRT = NBACK

                  END IF

               END IF

!              --- assign point IC to part that is currently being created

               IPOWN(IC) = LPARTS(1,JPART)
               NINPRT    = NINPRT + 1

            END IF

         END DO
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWORB ( IPOWN, IDIR, IPART, NPART, LPARTS,
     &                   MXC  , MYC )
!
!****************************************************************
!
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Performs an Orthogonal Recursive Bisection partitioning
!
!  3. Method
!
!     Starting with a single part (the entire domain), each part is
!     recursively partitioned by bisecting it, until all parts have been
!     created. The bisection direction is swapped in each direction.
!
!     This method is described in Ph.D. Thesis of M. Roest
!     entitled:
!     Partitioning for parallel finite difference computations
!     in coastal water simulation, DUT, 1997
!
!  4. Argument variables
!
!     IDIR        direction of cutting
!                 1 = row
!                 2 = column
!     IPART       part number that must be partitioned
!     IPOWN       array giving the subdomain number of each gridpoint
!     LPARTS      list of parts to be created
!                    lparts(1,i) = number of i-th part to be created
!                    lparts(2,i) = size of i-th part to be created
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!     NPART       number of parts to be created
!
      INTEGER   IDIR, IPART, MXC, MYC, NPART
      INTEGER   IPOWN(*)
      INTEGER*8 LPARTS(2,*)
!
!  6. Local variables
!
!     IDIFF :     the difference to be applied to a subdomain-size
!     IENT  :     number of entries
!     IP    :     counter of parts to be splitted
!     ISPLIT:     counter of parts to be created in splitting
!     ISSUCC:     flag indicating success in reducing a difference in size
!                 0=no
!                 1=yes
!     IWORK :     see description LPARTS
!     J     :     loop counter
!     JEND  :     number of last new part to be created by splitting
!     JPARTE:     number of last part in 1..npart belonging to jpart
!     JPARTS:     number of first part in 1..npart belonging to jpart
!     JSTART:     number of first new part to be created by splitting
!     KSPLIT:     number of parts to be created in a particular splitting
!     NP    :     number of parts to be created in a particular recursion
!     NSPLIT:     maximum number of parts to be created in one splitting
!                 (NB: 2 = bisection, 4 = quadrisection)

      INTEGER   IDIFF, IENT, IP, ISPLIT, ISSUCC, J, JEND, JPARTE,
     &          JPARTS, JSTART, KSPLIT, NP, NSPLIT
      INTEGER*8 IWORK(2,NPART)
      PARAMETER(NSPLIT = 2)
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     SWSTRIP          Performs a stripwise partitioning with straight
!                      interfaces
!
!  9. Subroutines calling
!
!     SWPARTIT
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
!     while not enough parts have been created, do another recursion
!
!        for each part that currently exists
!
!          determine which final parts belong to this part
!          if the number of such parts > 1, do further splitting
!
!            determine into how many parts this part must be split
!            determine the sizes and numbers of parts to be created
!
!            do splitting
!
!            determine whether objective partsizes have been modified
!            and distribute the difference over the constituent parts
!
!        swap cutting direction
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWORB')

!     --- while not enough parts have been created, do another recursion

      NP = 1
100   IF ( NP.LE.NPART ) THEN

!        --- for each part that currently exists

         DO IP = 1, NP

!           --- determine which final parts belong to this part

            JPARTS = (IP-1)*NPART/NP+1
            JPARTE = (IP  )*NPART/NP

!           --- if the number of such parts > 1, do further splitting

            IF ( (JPARTE-JPARTS+1).GT.1 ) THEN

!              --- determine into how many parts this part must be split

               KSPLIT = MIN(NSPLIT,JPARTE-JPARTS+1)

!              --- determine the sizes and numbers of parts to be created

               DO ISPLIT = 1, KSPLIT

                  JSTART = JPARTS+(ISPLIT-1)*(JPARTE-JPARTS+1)/KSPLIT
                  JEND   = JPARTS+    ISPLIT*(JPARTE-JPARTS+1)/KSPLIT-1

                  IWORK(1,ISPLIT) = JSTART

                  IWORK(2,ISPLIT) = 0
                  DO J = JSTART, JEND
                     IWORK(2,ISPLIT) = IWORK(2,ISPLIT) + LPARTS(2,J)
                  END DO

               END DO

!              --- do splitting

               CALL SWSTRIP (IPOWN,IDIR,JPARTS,KSPLIT,IWORK,MXC,MYC)

!              --- determine whether objective partsizes have been modified
!                  in SWSTRIP in order to make straight interfaces

               DO ISPLIT = 1, KSPLIT

                  JSTART = JPARTS+(ISPLIT-1)*(JPARTE-JPARTS+1)/KSPLIT
                  JEND   = JPARTS+    ISPLIT*(JPARTE-JPARTS+1)/KSPLIT-1

                  DO J = JSTART, JEND
                     IWORK(2,ISPLIT) = IWORK(2,ISPLIT) - LPARTS(2,J)
                  END DO

!                 --- and distribute the difference over the contiguous
!                     parts making sure not to cause negative subdomain-sizes

                  J      = JSTART
                  ISSUCC = 0
                  IF ( IWORK(2,ISPLIT).LT.0 ) THEN
                     IDIFF = -1
                  ELSE
                     IDIFF =  1
                  END IF

!                 --- reduce the difference until nothing is left

200               IF ( IWORK(2,ISPLIT).NE.0 ) THEN

!                    --- only adjust parts if their size remains valid

                     IF ( LPARTS(2,J).GT.0 .AND.
     &                    ((LPARTS(2,J)+IDIFF).GT.0) ) THEN
                        ISSUCC          = 1
                        LPARTS(2,J)     = LPARTS(2,J)     + IDIFF
                        IWORK(2,ISPLIT) = IWORK(2,ISPLIT) - IDIFF
                     END IF

!                    --- go on to the next part

                     J = J + 1

!                    --- when all parts have been visited, go back to first

                     IF ( J.GT.JEND ) THEN

!                       --- check whether any reduction of the difference
!                           was done in last pass over all parts

                        IF ( ISSUCC.EQ.0 ) THEN
                           CALL MSGERR (4,'Internal problem in SWORB')
                           RETURN
                        END IF
                        J = JSTART
                     END IF

                     GOTO 200
                  END IF

               END DO

            END IF

         END DO

!        --- swap cutting direction

         IDIR = MOD(IDIR,2) + 1

         NP = NSPLIT*NP
         GOTO 100

      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWPARTIT ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Sep. 04: determines load per processor based on speed
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Carries out the partitioning of the SWAN computational grid
!
!  3. Method
!
!     Based on Orthogonal Recursive Bisection
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  6. Local variables
!
!     I     :     loop counter
!     ICNT  :     auxiliary integer to count weights
!     IDIR  :     direction of cutting
!                 1 = row
!                 2 = column
!     IENT  :     number of entries
!     IX    :     index in x-direction
!     IY    :     index in y-direction
!     IWORK :     work array with the following meaning:
!                    IWORK(1,i) = number of i-th part to be created
!                    IWORK(2,i) = size of i-th part to be created
!     NACTP :     total number of active gridpoints
!     NPCUM :     cumulative number of gridpoints
!
      INTEGER   I, IDIR, IENT, IX, IY
      INTEGER*8 ICNT, NACTP, NPCUM
      INTEGER*8 IWORK(2,NPROC)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWORB            Performs an Orthogonal Recursive Bisection partitioning
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWDECOMP
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
!     if not parallel, return
!     determine direction of cutting
!     determine number of active points
!     determine numbers and sizes of parts to be created
!     partition grid
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWPARTIT')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- determine direction of cutting

      IF ( KPART.EQ.3 ) THEN
         IDIR = 1
      ELSE IF ( KPART.EQ.4 ) THEN
         IDIR = 2
      ELSE
         IF ( MXC.GT.MYC ) THEN
            IDIR = 2
            IF (.NOT.LORB) KPART = 4
         ELSE
            IDIR = 1
            IF (.NOT.LORB) KPART = 3
         END IF
      END IF

!     --- determine number of active points and
!         set IPOWN to 1 in these points

      NACTP = 0
      DO IX = 1, MXC
         DO IY = 1, MYC
            IF ( KGRPGL(IX,IY).NE.1 ) THEN
               IPOWN(IX,IY) = 1
               NACTP        = NACTP + 1
            END IF
         END DO
      END DO

!     --- determine numbers and sizes of parts to be created

      NPCUM = 0
      ICNT  = 0
      DO I = 1, NPROC
         ICNT       = ICNT + IWEIG(I)
         IWORK(1,I) = I
         IWORK(2,I) = (NACTP*ICNT)/SUM(IWEIG) - NPCUM
         NPCUM      = (NACTP*ICNT)/SUM(IWEIG)
      END DO
      DEALLOCATE(IWEIG)

!     --- partition grid

      IF ( LORB ) THEN

!        --- performs orthogonal recursive bisection partitioning

         CALL SWORB ( IPOWN, IDIR, 1, NPROC, IWORK, MXC, MYC )
         IF (STPNOW()) RETURN

      ELSE

!        --- performs stripwise partitioning

         CALL SWSTRIP ( IPOWN, IDIR, 1, NPROC, IWORK, MXC, MYC )
         IF (STPNOW()) RETURN

      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWBLADM ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Jul. 04: determine global bounds in subdomains
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     For the present node, carries out the block administration
!     and determines array bounds with respect to global grid
!
!  3. Method
!
!     Based on domain decomposition, the interface sizes are
!     determined that is needed for the setup of block
!     administration stored as IBLKAD
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     CHARS :     character for passing info to MSGERR
!     I     :     loop counter
!     ICNT  :     counter
!     ICRECV:     array containing positions of unknowns
!                 to be received from neighbour
!     ICSEND:     array containing positions of unknowns
!                 to be sent to neighbour
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IHALO :     actual width of halo area
!     IL    :     last non-character in string
!     INB   :     neighbour counter
!     ISTART:     startaddress for each size interface in array IBLKAD
!     ITEMP :     temporary work array to store block administration
!     IX    :     index in x-direction
!     IXOFF :     offset in x-direction
!     IY    :     index in y-direction
!     IYOFF :     offset in y-direction
!     IWORK :     array used to determine interface sizes
!                   IWORK(1,i) = number of the i-th neighbour
!                   IWORK(2,i) = position of the i-th neighbour with
!                                respect to present subdomain
!                                (resp. top, bottom, right, left)
!                   IWORK(3,i) = size of interface to i-th neighbour
!     J     :     loop counter
!     JE    :     last index of part of halo area
!     JS    :     first index of part of halo area
!     JOFFS :     offsets at which a point of a neigbhour domain can be found
!     LENGTH:     actual length of array IBLKAD
!     MSGSTR:     string to pass message to call MSGERR
!     MXSIZ :     size of present subdomain in x-direction
!     MYSIZ :     size of present subdomain in y-direction
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping points in halo area
!     NOVLUP:     number of overlapping points in part of halo area
!     NPCUM :     cumulative number of overlapping points
!
      INTEGER      I, ICNT, IDOM, IENT, IF, IHALO, IL, INB, ISTART,
     &             IX, IXOFF, IY, IYOFF, J, JE, JS, JOFFS(2,4), LENGTH,
     &             MXSIZ, MYSIZ, NNEIGH, NOVLU, NOVLUP, NPCUM
      INTEGER      IWORK(3,NPROC),
     &             ICRECV(NPROC,MAX(IHALOX,IHALOY)*MAX(MXC,MYC)),
     &             ICSEND(NPROC,MAX(IHALOX,IHALOY)*MAX(MXC,MYC))
      CHARACTER*20 INTSTR, CHARS
      CHARACTER*80 MSGSTR

      INTEGER, ALLOCATABLE :: ITEMP(:)
!
!  8. Subroutines used
!
!     INTSTR           Converts integer to string
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWDECOMP
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
!     intialize offsets to be used in searching for interfaces
!     determine enclosing box of present subdomain
!     if subdomain appears to be empty
!        give warning and set empty bounding box
!     else
!        extend enclosing box to include halo area
!     localize global bounds in present subdomain
!     determine size of enclosing box
!     determine interface sizes:
!
!        loop over global grid
!           if point belongs to this part
!              for each of the four sizes
!                  if a neighbouring subdomain is found there
!                     find it in the list of neighbours
!                     if not yet in the list, add it
!                     store position of neighbour
!                     update number of overlapping unknowns
!
!     store block administration
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWBLADM')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- intialize offsets to be used in searching for interfaces

      JOFFS = RESHAPE((/0,1,0,-1,1,0,-1,0/), (/2,4/))

!     --- determine enclosing box of present subdomain

      MXF = MXC+1
      MYF = MYC+1
      MXL = 0
      MYL = 0

      DO IX = 1, MXC
         DO IY = 1, MYC

            IF( IPOWN(IX,IY).EQ.INODE ) THEN

               MXF = MIN(IX,MXF)
               MYF = MIN(IY,MYF)
               MXL = MAX(IX,MXL)
               MYL = MAX(IY,MYL)

            END IF

         END DO
      END DO

!     --- if subdomain appears to be empty

      IF ( MXF.GT.MXL .OR. MYF.GT.MYL ) THEN

!        --- give warning and set empty bounding box

         CHARS = INTSTR(INODE)
         CALL TXPBLA(CHARS,IF,IL)
         MSGSTR = 'Empty subdomain is detected - '//
     &            ' node number is '//CHARS(IF:IL)
         CALL MSGERR ( 2, MSGSTR )

         MXF = 1
         MYF = 1
         MXL = 0
         MYL = 0

      ELSE

!        --- extend enclosing box to include halo area

         MXF = MAX(1  ,MXF-IHALOX)
         MYF = MAX(1  ,MYF-IHALOY)
         MXL = MIN(MXC,MXL+IHALOX)
         MYL = MIN(MYC,MYL+IHALOY)

      END IF

!     --- localize global bounds in present subdomain

      LMXF = MXF.EQ.1
      LMXL = MXL.EQ.MXCGL
      LMYF = MYF.EQ.1
      LMYL = MYL.EQ.MYCGL

!     --- determine size of enclosing box

      MXSIZ = MXL - MXF + 1
      MYSIZ = MYL - MYF + 1

      IWORK  = 0
      ICRECV = 0
      ICSEND = 0

!     --- determine interface sizes

      DO IX = 1, MXC
         DO IY = 1, MYC

!           --- if point belongs to this part

            IF ( IPOWN(IX,IY).EQ.INODE ) THEN

!              --- for each of the four sizes

               DO I = 1, 4

                  IXOFF = JOFFS(1,I)
                  IYOFF = JOFFS(2,I)

                  IF ( I.EQ.1 .OR. I.EQ.2 ) THEN
                     IHALO = IHALOY
                  ELSE
                     IHALO = IHALOX
                  END IF

!                 --- if a neighbouring subdomain is found there

                  IF ( (IX+IXOFF).GT.0.AND.(IX+IXOFF).LE.MXC.AND.
     &                 (IY+IYOFF).GT.0.AND.(IY+IYOFF).LE.MYC ) THEN

                     IF ( IPOWN(IX+IXOFF,IY+IYOFF).NE.0.AND.
     &                    IPOWN(IX+IXOFF,IY+IYOFF).NE.INODE ) THEN

!                       --- find it in the list of neighbours

                        IDOM = IPOWN(IX+IXOFF,IY+IYOFF)

                        INB = 1
  100                   IF ( INB.LE.NPROC .AND.
     &                       IWORK(1,INB).NE.IDOM .AND.
     &                       IWORK(1,INB).NE.0 )  THEN
                           INB = INB + 1
                           GOTO 100
                        END IF

                        IF ( INB.GT.NPROC ) THEN
                          CALL MSGERR (4,'Found more neighbours than '//
     &                                 'subdomains in the partitioning')
                          RETURN
                        END IF

!                       --- if not yet in the list, add it

                        IF ( IWORK(1,INB).EQ.0 ) IWORK(1,INB) = IDOM

!                       --- store position of neighbour with respect to
!                           present subdomain

                        IWORK(2,INB) = I

!                       --- update number of overlapping unknowns

                        IWORK(3,INB) = IWORK(3,INB) + IHALO

                        DO J = 1, IHALO
                           ICSEND(INB,IWORK(3,INB)-J+1) =
     &                                      (IY-(J-1)*IYOFF-MYF)*MXSIZ +
     &                                      (IX-(J-1)*IXOFF-MXF+1)
                           ICRECV(INB,IWORK(3,INB)-J+1) =
     &                                      (IY+    J*IYOFF-MYF)*MXSIZ +
     &                                      (IX+    J*IXOFF-MXF+1)
                        END DO

                     END IF

                  END IF

               END DO

            END IF

         END DO
      END DO

!     --- store block administration

      ALLOCATE(ITEMP(1+2*NPROC*(2+MAX(IHALOX,IHALOY)*MAX(MXC,MYC))))
      ITEMP = -999

      NNEIGH    = COUNT(IWORK(1,:)>0)
      ITEMP(1) = NNEIGH
      ISTART    = 3*NNEIGH+2
      DO INB = 1, NNEIGH
         ITEMP(3*INB-1) = IWORK(1,INB)
         ITEMP(3*INB  ) = IWORK(2,INB)
         ITEMP(3*INB+1) = ISTART
         NOVLU          = IWORK(3,INB)
         ITEMP(ISTART)  = NOVLU
         DO I = 1, NOVLU
            ITEMP(ISTART      +I) = ICSEND(INB,I)
            ITEMP(ISTART+NOVLU+I) = ICRECV(INB,I)
         END DO
         ISTART = ISTART + 2*NOVLU+1
      END DO

      LENGTH = COUNT(ITEMP/=-999)
      ALLOCATE(IBLKAD(LENGTH))

      DO I = 1, LENGTH
         IBLKAD(I) = ITEMP(I)
      END DO

      DEALLOCATE(ITEMP)

!     --- determine whether halo area must be partition
!         when dealing with large amount of data

      IF ( LORB ) THEN
         NHALOP = 1
      ELSE
         NOVLU = MAXVAL(IWORK(3,:))
         CALL SWREDUCE ( NOVLU, 1, SWINT, SWMAX )
         IF ( NOVLU.GT.MAXBUF ) THEN
            NHALOP = NOVLU / MAXBUF
         ELSE
            NHALOP = 1
         ENDIF
      END IF

      IF ( NHALOP.GT.10 ) THEN
         CALL MSGERR (1, 'halo area is divided in more than 10 parts')
      END IF

      IF ( NHALOP.GT.1 ) THEN
         ALLOCATE(JHALOE(NNEIGH,NHALOP))
         ALLOCATE(JHALOS(NNEIGH,NHALOP))
         DO INB = 1, NNEIGH
            NOVLU = IWORK(3,INB)
            NPCUM = 0
            ICNT  = 0
            JE    = 0
            DO I = 1, NHALOP
               ICNT          = ICNT + 1
               NOVLUP        = (NOVLU*ICNT)/NHALOP - NPCUM
               NPCUM         = (NOVLU*ICNT)/NHALOP
               JS            = JE + 1
               JE            = JS + NOVLUP - 1
               JHALOS(INB,I) = JS
               JHALOE(INB,I) = JE
            END DO
         END DO
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCORN ( IPOWN, MXC, MYC )
!
!****************************************************************
!
      USE OCPCOMM4
      USE M_PARALL
!
      IMPLICIT NONE
!
!
!  0. Authors
!
!      1.10: Marcel Zijlema
!
!  1. Updates
!
!      1.10, Feb. 12: New subroutine
!
!  2. Purpose
!
!     For own subdomain find the grid points in the
!     'halo' corners in the case of ORB partitioning
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     IPOWN       array giving the subdomain number of each gridpoint
!     MXC         maximum counter of gridpoints in x-direction
!     MYC         maximum counter of gridpoints in y-direction
!
      INTEGER MXC, MYC
      INTEGER IPOWN(MXC,MYC)
!
!  5. Parameter variables
!
!     ---
!
!  6. Local variables
!
!     I     :     loop counter
!     ICOR  :     index corner point in x-direction
!     IENT  :     number of entries
!     IX    :     index in x-direction
!     IXOFF :     offset in x-direction
!     IY    :     index in y-direction
!     IYOFF :     offset in y-direction
!     JCOR  :     index corner point in y-direction
!     JOFFS :     offsets at which a corner of subdomain can be found
!     MXSUB :     start and end indices of subdomain in x-direction
!     MYSUB :     start and end indices of subdomain in y-direction
!
      INTEGER      I, IENT, IX, IXOFF, IY, IYOFF
      INTEGER      ICOR(4,2), JCOR(4,2), JOFFS(2,4), MXSUB(2), MYSUB(2)
!
!  8. Subroutines used
!
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWDECOMP
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
!     find grid points in each 'halo' corner of own subdomain
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCORN')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- intialize offsets to be used in searching for corners

      JOFFS = RESHAPE((/1,1,1,-1,-1,-1,-1,1/),(/2,4/))

!     --- initialize some arrays

      ICOR(:,1) = MXF+1
      ICOR(:,2) = 0
      JCOR(:,1) = MYF+1
      JCOR(:,2) = 0

      MXSUB(1) = 1
      MXSUB(2) = MXCGL
      MYSUB(1) = 1
      MYSUB(2) = MYCGL
      IF (.NOT.LMXF) MXSUB(1) = MXF+IHALOX
      IF (.NOT.LMXL) MXSUB(2) = MXL-IHALOX
      IF (.NOT.LMYF) MYSUB(1) = MYF+IHALOY
      IF (.NOT.LMYL) MYSUB(2) = MYL-IHALOY

!     --- for each of the four corners

      DO I = 1, 4

         IXOFF = JOFFS(1,I)
         IYOFF = JOFFS(2,I)

         IX = MXSUB((IXOFF+3)/2)
         IY = MYSUB((IYOFF+3)/2)

!        --- if a corner is found there

         IF ( (IX+IXOFF).GT.0.AND.(IX+IXOFF).LE.MXC.AND.
     &        (IY+IYOFF).GT.0.AND.(IY+IYOFF).LE.MYC ) THEN

            IF ( IPOWN(IX+IXOFF,IY+IYOFF).NE.0.AND.
     &           IPOWN(IX+IXOFF,IY+IYOFF).NE.INODE ) THEN

!              --- determine start and end indices of the corner

               IF ( I.EQ.1 ) THEN
                  ICOR(1,1) = IX+IXOFF
                  ICOR(1,2) = MXL
                  JCOR(1,1) = IY+IYOFF
                  JCOR(1,2) = MYL
               ELSE IF ( I.EQ.2 ) THEN
                  ICOR(2,1) = IX+IXOFF
                  ICOR(2,2) = MXL
                  JCOR(2,1) = 1
                  JCOR(2,2) = IY+IYOFF
               ELSE IF ( I.EQ.3 ) THEN
                  ICOR(3,1) = 1
                  ICOR(3,2) = IX+IXOFF
                  JCOR(3,1) = 1
                  JCOR(3,2) = IY+IYOFF
               ELSE IF ( I.EQ.4 ) THEN
                  ICOR(4,1) = 1
                  ICOR(4,2) = IX+IXOFF
                  JCOR(4,1) = IY+IYOFF
                  JCOR(4,2) = MYL
               END IF

            END IF

         END IF

      END DO

!     --- determine for each grid point whether it is in the corner or not

      DO IX = MXF, MXL
         DO IY = MYF, MYL
            INCORN(IX,IY)= ( IX.GE.ICOR(1,1) .AND. IX.LE.ICOR(1,2) .AND.
     &                       IY.GE.JCOR(1,1) .AND. IY.LE.JCOR(1,2)) .OR.
     &                     ( IX.GE.ICOR(2,1) .AND. IX.LE.ICOR(2,2) .AND.
     &                       IY.GE.JCOR(2,1) .AND. IY.LE.JCOR(2,2)) .OR.
     &                     ( IX.GE.ICOR(3,1) .AND. IX.LE.ICOR(3,2) .AND.
     &                       IY.GE.JCOR(3,1) .AND. IY.LE.JCOR(3,2)) .OR.
     &                     ( IX.GE.ICOR(4,1) .AND. IX.LE.ICOR(4,2) .AND.
     &                       IY.GE.JCOR(4,1) .AND. IY.LE.JCOR(4,2))
         END DO
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWDECOMP
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata3                                                  40.41
      USE M_PARALL                                                        40.31
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
!      1.10: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!      1.10, Feb. 12: make it appropriate for ORB partitioning (for SWASH)
!
!  2. Purpose
!
!     Carries out domain decomposition meant for
!     distributed-memory approach
!
!  3. Method
!
!     First, carry out the partitioning of the
!     SWAN computational grid and then do the
!     block administration
!
!  4. Argument variables
!
!     ---
!
!  6. Local variables
!
!     IENT  :     number of entries
!     IPOWN :     array giving the subdomain number of each gridpoint
!
      INTEGER IENT
      INTEGER, ALLOCATABLE :: IPOWN(:,:)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWBLADM          Carries out the block administration
!     SWPARTIT         Carries out the partitioning of the SWAN
!                      computational grid
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWMAIN
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
!     store the original values of MXC, MYC and MCGRD
!     if not parallel, return
!     carry out the partitioning of computational grid
!     carry out the block administration
!     find corner points in case of ORB partitioning
!     compute MXC and MYC for each subdomain
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWDECOMP')

!     --- store the original values of MXC, MYC and MCGRD of
!         global computational grid

      MXCGL   = MXC
      MYCGL   = MYC
      MCGRDGL = MCGRD
      MXF     = 1
      MXL     = MXC
      MYF     = 1
      MYL     = MYC
      LMXF    = MXF.EQ.1                                                  40.41
      LMXL    = MXL.EQ.MXCGL                                              40.41
      LMYF    = MYF.EQ.1                                                  40.41
      LMYL    = MYL.EQ.MYCGL                                              40.41

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- carry out the partitioning of the SWAN
!         computational grid

      ALLOCATE(IPOWN(MXC,MYC))
      IPOWN = 0
      CALL SWPARTIT( IPOWN, MXC, MYC )
      IF (STPNOW()) RETURN

!     --- carry out the block administration

      CALL SWBLADM( IPOWN, MXC, MYC )
      IF (STPNOW()) RETURN

!     --- find corner points in case of ORB partitioning                   1.10
!         (they need to be excluded! See SwashInitCompGrid.f90)            1.10

      IF (LORB) CALL SWCORN( IPOWN, MXC, MYC )                             1.10

!     --- compute MXC and MYC for each subdomain

      MXC = MXL - MXF + 1
      MYC = MYL - MYF + 1

!     --- allocate work arrays for exchanging data

      ALLOCATE(IWRKEX(MAX(IHALOX,IHALOY)*MAX(MXC,MYC)))
      IWRKEX = 0
      ALLOCATE(WRKEXG(MAX(IHALOX,IHALOY)*MAX(MXC,MYC)))
      WRKEXG = 0.
      ALLOCATE(DWRKEX(MAX(IHALOX,IHALOY)*MAX(MXC,MYC)))
      DWRKEX = 0D0

      DEALLOCATE(IPOWN)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXCHG ( FIELD, KGRPNT, KS, KE )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata3                                                  40.41
      USE M_PARALL                                                        40.31
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
!     41.50: Marcel Zijlema
!      1.19: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.50, Aug. 10: extended with vertical direction for SWASH
!      1.19, Nov. 12: make exchanging large amount of data possible
!
!  2. Purpose
!
!     Updates field array through exchanging
!     values between neighbouring subdomains
!
!  3. Method
!
!     Made use of MPI by means of SWSENDNB and SWRECVNB
!     and also block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     FIELD       field array for which 'halo' values must
!                 be copied from neighbouring subdomains
!     KGRPNT      indirect addressing for grid points
!     KE          last index in vertical direction
!     KS          first index in vertical direction
!
      INTEGER KS, KE
      INTEGER KGRPNT(MXC*MYC)
      REAL    FIELD(MCGRD,KS:KE)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     ISTART:     pointer in array IBLKAD
!     ITAG  :     message tag for sending and receiving
!     J     :     loop counter
!     JE    :     last index of J
!     JS    :     first index of J
!     K     :     loop counter in vertical direction
!     N     :     loop counter
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!     NOVLUP:     number of overlapping unknowns
!                 in part of halo area
!
      INTEGER IDOM, IENT, INB, ISTART, ITAG, J, JE, JS, K,
     &        N, NNEIGH, NOVLU, NOVLUP
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!     SWSENDNB         Data is sent to a neighbour
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!
      LOGICAL STPNOW
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
!     if not parallel, return
!
!     consider either whole halo area or parts of it, do
!
!        for all layers do
!
!           for all neighbouring subdomains do
!              get subdomain number, pointer and size
!              store data to be sent in array WRKEXG
!              send array WRKEXG
!
!           for all neighbouring subdomains do
!              get subdomain number, pointer and size
!              receive next array and store in WRKEXG
!              store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWEXCHG')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(203)

      NNEIGH = IBLKAD(1)

      IF ( NHALOP.EQ.1 ) THEN
!        take halo area at once for exchanging data

!        --- for all layers do

         DO K = KS, KE

!           --- for all neighbouring subdomains do

            DO INB = 1, NNEIGH

!              --- get subdomain number, pointer and size

               IDOM   = IBLKAD(3*INB-1)
               ISTART = IBLKAD(3*INB+1)
               NOVLU  = IBLKAD(ISTART)

!              --- store data to be sent in array WRKEXG

               DO J = 1, NOVLU
                  WRKEXG(J) = FIELD(KGRPNT(IBLKAD(ISTART+J)),K)
               END DO

!              --- send array WRKEXG

               ITAG = 2
               CALL SWSENDNB ( WRKEXG, NOVLU, SWREAL, IDOM, ITAG )
               IF (STPNOW()) RETURN

            END DO

!           --- for all neighbouring subdomains do

            DO INB = 1, NNEIGH

!              --- get subdomain number, pointer and size

               IDOM   = IBLKAD(3*INB-1)
               ISTART = IBLKAD(3*INB+1)
               NOVLU  = IBLKAD(ISTART)

!              --- receive next array and store in WRKEXG

               ITAG = 2
               CALL SWRECVNB ( WRKEXG, NOVLU, SWREAL, IDOM, ITAG )
               IF (STPNOW()) RETURN

!              --- store the received data

               DO J = 1, NOVLU
                  FIELD(KGRPNT(IBLKAD(ISTART+NOVLU+J)),K) = WRKEXG(J)
               END DO

            END DO

         END DO

      ELSE
!        consider each part of halo area for exchanging data

!        --- for all layers do

         DO K = KS, KE

!           --- for all parts of halo area

            DO N = 1, NHALOP

!              --- for all neighbouring subdomains do

               DO INB = 1, NNEIGH

!                 --- get subdomain number, pointer and size

                  IDOM   = IBLKAD(3*INB-1)
                  ISTART = IBLKAD(3*INB+1)

                  JS = JHALOS(INB,N)
                  JE = JHALOE(INB,N)

                  NOVLUP = JE - JS + 1

!                 --- store data to be sent in array WRKEXG

                  DO J = JS, JE
                     WRKEXG(J-JS+1) = FIELD(KGRPNT(IBLKAD(ISTART+J)),K)
                  END DO

!                 --- send array WRKEXG

                  ITAG = 2
                  CALL SWSENDNB ( WRKEXG, NOVLUP, SWREAL, IDOM, ITAG )
                  IF (STPNOW()) RETURN

               END DO

!              --- for all neighbouring subdomains do

               DO INB = 1, NNEIGH

!                 --- get subdomain number, pointer and size

                  IDOM   = IBLKAD(3*INB-1)
                  ISTART = IBLKAD(3*INB+1)
                  NOVLU  = IBLKAD(ISTART)

                  JS = JHALOS(INB,N)
                  JE = JHALOE(INB,N)

                  NOVLUP = JE - JS + 1

!                 --- receive next array and store in WRKEXG

                  ITAG = 2
                  CALL SWRECVNB ( WRKEXG, NOVLUP, SWREAL, IDOM, ITAG )
                  IF (STPNOW()) RETURN

!                 --- store the received data

                  DO J = JS, JE
                     FIELD(KGRPNT(IBLKAD(ISTART+NOVLU+J)),K) =
     &                                                    WRKEXG(J-JS+1)
                  END DO

               END DO

            END DO

         END DO

      END IF

!TIMG      CALL SWTSTO(203)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXCHGI ( IFLD, KGRPNT, KS, KE )
!
!****************************************************************
!
      USE OCPCOMM4
      USE SwashCommdata3
      USE M_PARALL
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
!      1.00: Marcel Zijlema
!
!  1. Updates
!
!      1.00, Aug. 14: New subroutine
!
!  2. Purpose
!
!     Updates field array of type integer through exchanging
!     values between neighbouring subdomains
!
!  3. Method
!
!     Made use of MPI by means of SWSENDNB and SWRECVNB
!     and also block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     IFLD        field array for which 'halo' values must
!                 be copied from neighbouring subdomains
!     KGRPNT      indirect addressing for grid points
!     KE          last index in vertical direction
!     KS          first index in vertical direction
!
      INTEGER KS, KE
      INTEGER KGRPNT(MXC*MYC)
      INTEGER IFLD(MCGRD,KS:KE)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     ISTART:     pointer in array IBLKAD
!     ITAG  :     message tag for sending and receiving
!     J     :     loop counter
!     JE    :     last index of J
!     JS    :     first index of J
!     K     :     loop counter in vertical direction
!     N     :     loop counter
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!     NOVLUP:     number of overlapping unknowns
!                 in part of halo area
!
      INTEGER IDOM, IENT, INB, ISTART, ITAG, J, JE, JS, K,
     &        N, NNEIGH, NOVLU, NOVLUP
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!     SWSENDNB         Data is sent to a neighbour
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!
      LOGICAL STPNOW
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
!     if not parallel, return
!
!     consider either whole halo area or parts of it, do
!
!        for all layers do
!
!           for all neighbouring subdomains do
!              get subdomain number, pointer and size
!              store data to be sent in array IWRKEX
!              send array IWRKEX
!
!           for all neighbouring subdomains do
!              get subdomain number, pointer and size
!              receive next array and store in IWRKEX
!              store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWEXCHGI')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(203)

      NNEIGH = IBLKAD(1)

      IF ( NHALOP.EQ.1 ) THEN
!        take halo area at once for exchanging data

!        --- for all layers do

         DO K = KS, KE

!           --- for all neighbouring subdomains do

            DO INB = 1, NNEIGH

!              --- get subdomain number, pointer and size

               IDOM   = IBLKAD(3*INB-1)
               ISTART = IBLKAD(3*INB+1)
               NOVLU  = IBLKAD(ISTART)

!              --- store data to be sent in array IWRKEX

               DO J = 1, NOVLU
                  IWRKEX(J) = IFLD(KGRPNT(IBLKAD(ISTART+J)),K)
               END DO

!              --- send array IWRKEX

               ITAG = 2
               CALL SWSENDNB ( IWRKEX, NOVLU, SWINT, IDOM, ITAG )
               IF (STPNOW()) RETURN

            END DO

!           --- for all neighbouring subdomains do

            DO INB = 1, NNEIGH

!              --- get subdomain number, pointer and size

               IDOM   = IBLKAD(3*INB-1)
               ISTART = IBLKAD(3*INB+1)
               NOVLU  = IBLKAD(ISTART)

!              --- receive next array and store in IWRKEX

               ITAG = 2
               CALL SWRECVNB ( IWRKEX, NOVLU, SWINT, IDOM, ITAG )
               IF (STPNOW()) RETURN

!              --- store the received data

               DO J = 1, NOVLU
                  IFLD(KGRPNT(IBLKAD(ISTART+NOVLU+J)),K) = IWRKEX(J)
               END DO

            END DO

         END DO

      ELSE
!        consider each part of halo area for exchanging data

!        --- for all layers do

         DO K = KS, KE

!           --- for all parts of halo area

            DO N = 1, NHALOP

!              --- for all neighbouring subdomains do

               DO INB = 1, NNEIGH

!                 --- get subdomain number, pointer and size

                  IDOM   = IBLKAD(3*INB-1)
                  ISTART = IBLKAD(3*INB+1)

                  JS = JHALOS(INB,N)
                  JE = JHALOE(INB,N)

                  NOVLUP = JE - JS + 1

!                 --- store data to be sent in array IWRKEX

                  DO J = JS, JE
                     IWRKEX(J-JS+1) = IFLD(KGRPNT(IBLKAD(ISTART+J)),K)
                  END DO

!                 --- send array IWRKEX

                  ITAG = 2
                  CALL SWSENDNB ( IWRKEX, NOVLUP, SWINT, IDOM, ITAG )
                  IF (STPNOW()) RETURN

               END DO

!              --- for all neighbouring subdomains do

               DO INB = 1, NNEIGH

!                 --- get subdomain number, pointer and size

                  IDOM   = IBLKAD(3*INB-1)
                  ISTART = IBLKAD(3*INB+1)
                  NOVLU  = IBLKAD(ISTART)

                  JS = JHALOS(INB,N)
                  JE = JHALOE(INB,N)

                  NOVLUP = JE - JS + 1

!                 --- receive next array and store in IWRKEX

                  ITAG = 2
                  CALL SWRECVNB ( IWRKEX, NOVLUP, SWINT, IDOM, ITAG )
                  IF (STPNOW()) RETURN

!                 --- store the received data

                  DO J = JS, JE
                     IFLD(KGRPNT(IBLKAD(ISTART+NOVLU+J)),K) =
     &                                                    IWRKEX(J-JS+1)
                  END DO

               END DO

            END DO

         END DO

      END IF

!TIMG      CALL SWTSTO(203)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWEXCHGD ( FIELD, KGRPNT, KS, KE )
!
!****************************************************************
!
      USE OCPCOMM4
      USE SwashCommdata3
      USE M_PARALL
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
!      1.00: Marcel Zijlema
!
!  1. Updates
!
!      1.00, Aug. 14: New subroutine
!
!  2. Purpose
!
!     Updates field array of type double precision through exchanging
!     values between neighbouring subdomains
!
!  3. Method
!
!     Made use of MPI by means of SWSENDNBD and SWRECVNBD
!     and also block administration (stored in IBLKAD)
!
!  4. Argument variables
!
!     FIELD       field array for which 'halo' values must
!                 be copied from neighbouring subdomains
!     KGRPNT      indirect addressing for grid points
!     KE          last index in vertical direction
!     KS          first index in vertical direction
!
      INTEGER KS, KE
      INTEGER KGRPNT(MXC*MYC)
      REAL*8  FIELD(MCGRD,KS:KE)
!
!  6. Local variables
!
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     ISTART:     pointer in array IBLKAD
!     ITAG  :     message tag for sending and receiving
!     J     :     loop counter
!     JE    :     last index of J
!     JS    :     first index of J
!     K     :     loop counter in vertical direction
!     N     :     loop counter
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!     NOVLUP:     number of overlapping unknowns
!                 in part of halo area
!
      INTEGER IDOM, IENT, INB, ISTART, ITAG, J, JE, JS, K,
     &        N, NNEIGH, NOVLU, NOVLUP
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNBD        Data is received from a neighbour (double precision)
!     SWSENDNBD        Data is sent to a neighbour (double precision)
!TIMG!     SWTSTA           Start timing for a section of code
!TIMG!     SWTSTO           Stop timing for a section of code
!
      LOGICAL STPNOW
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
!     if not parallel, return
!
!     consider either whole halo area or parts of it, do
!
!        for all layers do
!
!           for all neighbouring subdomains do
!              get subdomain number, pointer and size
!              store data to be sent in array DWRKEX
!              send array DWRKEX
!
!           for all neighbouring subdomains do
!              get subdomain number, pointer and size
!              receive next array and store in DWRKEX
!              store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWEXCHGD')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!TIMG      CALL SWTSTA(203)

      NNEIGH = IBLKAD(1)

      IF ( NHALOP.EQ.1 ) THEN
!        take halo area at once for exchanging data

!        --- for all layers do

         DO K = KS, KE

!           --- for all neighbouring subdomains do

            DO INB = 1, NNEIGH

!              --- get subdomain number, pointer and size

               IDOM   = IBLKAD(3*INB-1)
               ISTART = IBLKAD(3*INB+1)
               NOVLU  = IBLKAD(ISTART)

!              --- store data to be sent in array DWRKEX

               DO J = 1, NOVLU
                  DWRKEX(J) = FIELD(KGRPNT(IBLKAD(ISTART+J)),K)
               END DO

!              --- send array DWRKEX

               ITAG = 2
               CALL SWSENDNBD ( DWRKEX, NOVLU, IDOM, ITAG )
               IF (STPNOW()) RETURN

            END DO

!           --- for all neighbouring subdomains do

            DO INB = 1, NNEIGH

!              --- get subdomain number, pointer and size

               IDOM   = IBLKAD(3*INB-1)
               ISTART = IBLKAD(3*INB+1)
               NOVLU  = IBLKAD(ISTART)

!              --- receive next array and store in DWRKEX

               ITAG = 2
               CALL SWRECVNBD ( DWRKEX, NOVLU, IDOM, ITAG )
               IF (STPNOW()) RETURN

!              --- store the received data

               DO J = 1, NOVLU
                  FIELD(KGRPNT(IBLKAD(ISTART+NOVLU+J)),K) = DWRKEX(J)
               END DO

            END DO

         END DO

      ELSE
!        consider each part of halo area for exchanging data

!        --- for all layers do

         DO K = KS, KE

!           --- for all parts of halo area

            DO N = 1, NHALOP

!              --- for all neighbouring subdomains do

               DO INB = 1, NNEIGH

!                 --- get subdomain number, pointer and size

                  IDOM   = IBLKAD(3*INB-1)
                  ISTART = IBLKAD(3*INB+1)

                  JS = JHALOS(INB,N)
                  JE = JHALOE(INB,N)

                  NOVLUP = JE - JS + 1

!                 --- store data to be sent in array DWRKEX

                  DO J = JS, JE
                     DWRKEX(J-JS+1) = FIELD(KGRPNT(IBLKAD(ISTART+J)),K)
                  END DO

!                 --- send array DWRKEX

                  ITAG = 2
                  CALL SWSENDNBD ( DWRKEX, NOVLUP, IDOM, ITAG )
                  IF (STPNOW()) RETURN

               END DO

!              --- for all neighbouring subdomains do

               DO INB = 1, NNEIGH

!                 --- get subdomain number, pointer and size

                  IDOM   = IBLKAD(3*INB-1)
                  ISTART = IBLKAD(3*INB+1)
                  NOVLU  = IBLKAD(ISTART)

                  JS = JHALOS(INB,N)
                  JE = JHALOE(INB,N)

                  NOVLUP = JE - JS + 1

!                 --- receive next array and store in DWRKEX

                  ITAG = 2
                  CALL SWRECVNBD ( DWRKEX, NOVLUP, IDOM, ITAG )
                  IF (STPNOW()) RETURN

!                 --- store the received data

                  DO J = JS, JE
                     FIELD(KGRPNT(IBLKAD(ISTART+NOVLU+J)),K) =
     &                                                    DWRKEX(J-JS+1)
                  END DO

               END DO

            END DO

         END DO

      END IF

!TIMG      CALL SWTSTO(203)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWRECVFLD ( FIELD, NFLD, IS, J, SWPDIR, KGRPNT )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata3                                                  40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Receives field arrays from bottom/top or left/right neighbour
!
!  3. Method
!
!     Use of SWRECVNB
!
!  4. Argument variables
!
!     FIELD       field arrays
!     IS          start index of J-th row
!     J           J-th row
!     KGRPNT      indirect addressing for grid points
!     NFLD        number of field arrays
!     SWPDIR      sweep direction
!
      INTEGER IS, J, NFLD, SWPDIR
      INTEGER KGRPNT(MXC,MYC)
      REAL    FIELD(MCGRD,NFLD)
!
!  6. Local variables
!
!     I     :     counter
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     IPNB  :     position of neighbour (=top, bottom, right, left)
!     ITAG  :     message tag for sending and receiving
!     K     :     loop counter
!     N     :     loop counter
!     NBOFF :     offset in neighbour position depending on sweep direction
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!     S     :     sign
!     WORK  :     work array to store data to received from neighbour
!
      INTEGER I, IDOM, IENT, INB, IPNB, ITAG, K, N, NBOFF, NNEIGH,
     &        NOVLU, S
      REAL    WORK(NFLD*MAX(IHALOX,IHALOY))
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWRECVNB         Data is received from a neighbour
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP
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
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get position
!        if neighbour is bottom/top or left/right
!           get subdomain number
!           receive next array and store in WORK
!           store the received data
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWRECVFLD')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IF ( SWPDIR.EQ.1 ) THEN
!        --- forward sweep
         NBOFF =  0
         S     = -1
      ELSE IF ( SWPDIR.EQ.2 ) THEN
!        --- backward sweep
         NBOFF = -1
         S     = +1
      END IF

      NOVLU = MAX(IHALOX,IHALOY)

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get position

         IPNB = IBLKAD(3*INB)

!        --- if neighbour is bottom/top or left/right

         IF ( IPNB.EQ.2+NBOFF .OR. IPNB.EQ.4+NBOFF ) THEN

!           --- get subdomain number

            IDOM = IBLKAD(3*INB-1)

!           --- receive next array and store in WORK

            ITAG  = 2
            CALL SWRECVNB ( WORK, NFLD*NOVLU, SWREAL, IDOM, ITAG )
            IF (STPNOW()) RETURN

!           --- store the received data

            I = 0
            IF ( KPART.EQ.4 ) THEN
               DO N = 1, NFLD
                  DO K = NOVLU, 1, -1
                     I = I + 1
                     FIELD(KGRPNT(IS+S*K,J),N) = WORK(I)
                  END DO
               END DO
            ELSE
               DO N = 1, NFLD
                  DO K = NOVLU, 1, -1
                     I = I + 1
                     FIELD(KGRPNT(J,IS+S*K),N) = WORK(I)
                  END DO
               END DO
            END IF

         END IF

      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWSENDFLD ( FIELD, NFLD, IE, J, SWPDIR, KGRPNT )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata3                                                  40.41
      USE M_PARALL                                                        40.31
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
!
!  1. Updates
!
!     40.30, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Sends field arrays to top/bottom or right/left neighbour
!
!  3. Method
!
!     Use of SWSENDNB
!
!  4. Argument variables
!
!     FIELD       field arrays
!     IE          end index of J-th row
!     J           J-th row
!     KGRPNT      indirect addressing for grid points
!     NFLD        number of field arrays
!     SWPDIR      sweep direction
!
      INTEGER IE, J, NFLD, SWPDIR
      INTEGER KGRPNT(MXC,MYC)
      REAL    FIELD(MCGRD,NFLD)
!
!  6. Local variables
!
!     I     :     counter
!     IDOM  :     subdomain number
!     IENT  :     number of entries
!     INB   :     neighbour counter
!     IPNB  :     position of neighbour (=top, bottom, right, left)
!     ITAG  :     message tag for sending and receiving
!     K     :     loop counter
!     N     :     loop counter
!     NBOFF :     offset in neighbour position depending on sweep direction
!     NNEIGH:     number of neighbouring subdomains
!     NOVLU :     number of overlapping unknowns
!     S     :     sign
!     WORK  :     work array to store data to be sent to neighbour
!
      INTEGER I, IDOM, IENT, INB, IPNB, ITAG, K, N, NBOFF, NNEIGH,
     &        NOVLU, S
      REAL    WORK(NFLD*MAX(IHALOX,IHALOY))
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWSENDNB         Data is sent to a neighbour
!
      LOGICAL STPNOW
!
!  9. Subroutines calling
!
!     SWCOMP
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
!     if not parallel, return
!
!     for all neighbouring subdomains do
!        get position
!        if neighbour is top/bottom or right/left
!           get subdomain number
!           store data to be sent in array WORK
!           send array WORK
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWSENDFLD')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

      IF ( SWPDIR.EQ.1 ) THEN
!        --- forward sweep
         NBOFF =  0
         S     = -1
      ELSE IF ( SWPDIR.EQ.2 ) THEN
!        --- backward sweep
         NBOFF =  1
         S     = +1
      END IF

      NOVLU = MAX(IHALOX,IHALOY)

      NNEIGH = IBLKAD(1)

!     --- for all neighbouring subdomains do

      DO INB = 1, NNEIGH

!        --- get position

         IPNB = IBLKAD(3*INB)

!        --- if neighbour is top/bottom or right/left

         IF ( IPNB.EQ.1+NBOFF .OR. IPNB.EQ.3+NBOFF ) THEN

!           --- get subdomain number

            IDOM = IBLKAD(3*INB-1)

!           --- store data to be sent in array WORK

            I = 0
            IF ( KPART.EQ.4 ) THEN
               DO N = 1, NFLD
                  DO K = NOVLU-1, 0, -1
                     I = I + 1
                     WORK(I) = FIELD(KGRPNT(IE+S*K,J),N)
                  END DO
               END DO
            ELSE
               DO N = 1, NFLD
                  DO K = NOVLU-1, 0, -1
                     I = I + 1
                     WORK(I) = FIELD(KGRPNT(J,IE+S*K),N)
                  END DO
               END DO
            END IF

!           --- send array WORK

            ITAG = 2
            CALL SWSENDNB ( WORK, NFLD*NOVLU, SWREAL, IDOM, ITAG )
            IF (STPNOW()) RETURN

         END IF

      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLLECT ( FIELDGL, FIELD )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata3                                                  40.41
      USE m_genarr                                                        40.31
      USE M_PARALL                                                        40.31
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
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.30, Mar. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.51, Feb. 05: extended to full arrays
!
!  2. Purpose
!
!     Collects FIELD array from all nodes
!
!  3. Method
!
!     Made use of MPI by means of SWGATHER
!
!  4. Argument variables
!
!     FIELD       field array in own subdomain
!     FIELDGL     global field array gathered from all nodes
!
      REAL    FIELD(*), FIELDGL(*)
!
!  6. Local variables
!
!     FLDC  :     auxiliary array for collecting data
!     IARRC :     auxiliary array for collecting grid indices
!     IARRL :     auxiliary array containing grid indices of each subdomain
!     IENT  :     number of entries
!     ILEN  :     length of field of current subdomain
!     ILENC :     length of field containing collected data
!     INDX  :     pointer in array
!     INDXC :     pointer in collected array
!     IP    :     node number
!     ISTART:     start pointer for each subdomain in array FLDC
!     IX    :     loop counter
!     IY    :     loop counter
!     MXFGL :     first index w.r.t. global grid in x-direction
!     MXLGL :     last index w.r.t. global grid in x-direction
!     MXSIZ :     size of present subdomain in x-direction
!     MYFGL :     first index w.r.t. global grid in y-direction
!     MYLGL :     last index w.r.t. global grid in y-direction
!
      INTEGER IENT, ILEN, ILENC, INDX, INDXC, IP, ISTART, IX, IY, MXSIZ,
     &        MXFGL, MXLGL, MYFGL, MYLGL
      INTEGER IARRL(4), IARRC(4,0:NPROC-1)

      REAL,    ALLOCATABLE :: FLDC(:)
!
!  8. Subroutines used
!
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWGATHER         Gathers different amounts of data from all nodes
!
      LOGICAL STPNOW
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
!     if not parallel, return
!     gather grid indices of all subdomains
!     determine total length for collecting data of all nodes and allocate array
!     gather data of all subdomains
!     copy gathered data to global array in appropriate manner
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLLECT')

!     --- if not parallel, return
      IF (.NOT.PARLL) RETURN

!     --- gather grid indices of all subdomains

      IARRL(1) = MXF
      IARRL(2) = MXL
      IARRL(3) = MYF
      IARRL(4) = MYL
      CALL SWGATHER (IARRC, 4*NPROC, IARRL, 4, SWINT )
      IF (STPNOW()) RETURN

!     --- determine total length for collecting data of all nodes and allocate array

      IF (INODE.EQ.MASTER) THEN
         ILENC = 0
         DO IP = 0, NPROC-1
            ILENC = ILENC +
     &           (IARRC(2,IP)-IARRC(1,IP)+1)*(IARRC(4,IP)-IARRC(3,IP)+1)
         END DO
         ALLOCATE(FLDC(ILENC))
      END IF
      ILEN = MXC*MYC

!     --- gather data of all subdomains

      CALL SWGATHER ( FLDC, ILENC, FIELD, ILEN, SWREAL )
      IF (STPNOW()) RETURN

!     --- copy gathered data to global array in appropriate manner

      IF (INODE.EQ.MASTER) THEN

         ISTART = 0

         DO IP = 0, NPROC-1
            IF ( IARRC(1,IP).EQ.1 ) THEN
               MXFGL = 1
            ELSE
               MXFGL = IARRC(1,IP) + IHALOX
            END IF
            IF ( IARRC(2,IP).EQ.MXCGL ) THEN
               MXLGL = MXCGL
            ELSE
               MXLGL = IARRC(2,IP) - IHALOX
            END IF
            IF ( IARRC(3,IP).EQ.1 ) THEN
               MYFGL = 1
            ELSE
               MYFGL = IARRC(3,IP) + IHALOY
            END IF
            IF ( IARRC(4,IP).EQ.MYCGL ) THEN
               MYLGL = MYCGL
            ELSE
               MYLGL = IARRC(4,IP) - IHALOY
            END IF

            MXSIZ = IARRC(2,IP)-IARRC(1,IP)+1

            DO IX = MXFGL, MXLGL
               DO IY = MYFGL, MYLGL
                  INDX  = (IY-1)*MXCGL+IX
                  INDXC = ISTART + (IY-IARRC(3,IP))*MXSIZ +
     &                              IX-IARRC(1,IP)+1
                  FIELDGL(INDX)= FLDC(INDXC)
               END DO
            END DO

            ISTART = ISTART + (IARRC(2,IP)-IARRC(1,IP)+1)*
     &                        (IARRC(4,IP)-IARRC(3,IP)+1)

         END DO

      END IF

      IF (INODE.EQ.MASTER) DEALLOCATE(FLDC)

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLOUT ( BLKND )                                       40.51
!
!****************************************************************
!
      USE SwashTimecomm                                                   40.41
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE SwashCommdata2                                                  40.51
      USE SwashCommdata3                                                  40.41
      USE outp_data
      USE M_PARALL                                                        40.31
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
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: re-design output process in parallel mode
!
!  2. Purpose
!
!     Collects output results
!
!  3. Method
!
!     Read individual process output files containing tables,
!     and block data and write them to
!     generic output files in appropriate manner
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain        40.41
!
      REAL    BLKND(MXCGL,MYCGL)                                          40.51
!
!  6. Local variables
!
!     CORQ  :     current item in list of request outputs
!     CUOPS :     current item in list of point sets
!     DTTIWR:     to write time string
!     IC    :     loop variable
!     IENT  :     number of entries
!     IP    :     loop variable
!     IRQ   :     request number
!     IT    :     time step counter
!     IT0   :     integer indicating first step of simulation
!     IT1   :     integer indicating last step of simulation
!     ITMP1 :     auxiliary integer
!     ITMP2 :     auxiliary integer
!     ITMP3 :     auxiliary integer
!     ITMP4 :     auxiliary integer
!     ITMP5 :     auxiliary integer
!     ITMP6 :     auxiliary integer
!     IUNIT :     counter for file unit numbers
!     MIP   :     total number of output points
!     MXK   :     number of points in x-direction of output frame
!     MYK   :     number of points in y-direction of output frame
!     OPENED:     logical whether a file is open or not
!     PSTYPE:     type of point set
!     RTMP1 :     auxiliary real
!     RTMP2 :     auxiliary real
!     RTMP3 :     auxiliary real
!     RTMP4 :     auxiliary real
!     RTYPE :     type of request
!     SNAMPF:     name of plot frame
!     TNEXT :     time of next requested output
!     XC    :     computational grid x-coordinate of output point
!     XP    :     user x-coordinate of output point
!     YC    :     computational grid y-coordinate of output point
!     YP    :     user y-coordinate of output point
!
      INTEGER   IC, IENT, IP, IRQ, IT, IT0, IT1, IUNIT, MIP, MXK, MYK
      INTEGER   ITMP1, ITMP2, ITMP3, ITMP4, ITMP5, ITMP6
      REAL*8    TNEXT
      REAL      RTMP1, RTMP2, RTMP3, RTMP4
      REAL, ALLOCATABLE :: XC(:), YC(:), XP(:), YP(:)                     40.51
      LOGICAL   OPENED
      CHARACTER PSTYPE*1, RTYPE*4, SNAMPF*8, DTTIWR*18
      TYPE(ORQDAT), POINTER :: CORQ
      TYPE(OPSDAT), POINTER :: CUOPS
!
!  8. Subroutines used
!
!     EQREAL           Logical comparing two reals
!     STPNOW           Logical indicating whether program must
!                      terminated or not
!     STRACE           Tracing routine for debugging
!     SWCOLBLK         Collects block output
!     SWCOLTAB         Collects table ouput
!     SwashCoorOutp    Computes coordinates of output points              40.51
!
      LOGICAL   EQREAL, STPNOW
!
!  9. Subroutines calling
!
!     SWMAIN
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
!     do for all COMPUTE commands
!        do for all time steps
!           do for all output requests
!              processing of output instructions necessary for collection
!              check time of output action
!              compute coordinates of output points
!              correct problem coordinates with offset values
!              rewrite table output by means of collection of output
!              rewrite block output by means of collection of output
!     close all files and delete process files
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLOUT')

      IF ( NREOQ.EQ.0 ) RETURN

!     --- do for all COMPUTE commands

      DO IC = 1, NCOMPT

         NSTATC = NINT(RCOMPT(IC,1))
         IF ( NSTATC.EQ.1 ) THEN
            IT0 = 0
         ELSE
            IT0 = 1
         END IF
         TFINC = RCOMPT(IC,3)
         TINIC = RCOMPT(IC,4)
         DT    = ETORQ(2)
         IT1   = NINT((TFINC-TINIC)/DT)
         TIMCO = TINIC

!        --- do for all time steps

         DO IT = IT0, IT1

            IF (NSTATM.GT.0) CHTIME = DTTIWR(ITMOPT, TIMCO)

!           --- do for all output requests

            CORQ => FORQ
            DO 100 IRQ = 1, NREOQ

!              --- processing of output instructions necessary
!                  for collection

!              --- check time of output action

               IF ( IT.EQ.IT0 .AND. IC.EQ.1 ) THEN
                  CORQ%OQR(1) = OURQT(IRQ)                                40.51
               END IF
               IF (CORQ%OQR(1).LT.TINIC) THEN
                  TNEXT = TINIC
               ELSE
                  TNEXT = CORQ%OQR(1)
               ENDIF
               IF ( REAL(TNEXT).GT.REAL(ETORQ(1)) ) THEN
                  GOTO 50
               ELSEIF (IT.EQ.IT1 .AND. CORQ%OQR(2).LT.0.) THEN
                  CORQ%OQR(1) = TIMCO
               ELSE IF ( CORQ%OQR(2).GT.0. .AND.
     &                  .NOT.REAL(TIMCO).LT.REAL(TNEXT) ) THEN
                  CORQ%OQR(1) = TNEXT + CORQ%OQR(2)
               ELSE
                  GOTO 50
               END IF

               RTYPE  = CORQ%RQTYPE
               SNAMPF = CORQ%PSNAME

               CUOPS => FOPS
               DO
                 IF (CUOPS%PSNAME.EQ.SNAMPF) EXIT
                 IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) GOTO 50
                 CUOPS => CUOPS%NEXTOPS
               END DO
               PSTYPE = CUOPS%PSTYPE

               IF ( PSTYPE.EQ.'F' .OR. PSTYPE.EQ.'H' ) THEN
                  MXK = CUOPS%OPI(1)
                  MYK = CUOPS%OPI(2)
                  MIP = MXK * MYK
               ELSE IF ( PSTYPE.EQ.'C' .OR. PSTYPE.EQ.'P' .OR.
     &                   PSTYPE.EQ.'N' ) THEN
                  MXK = 0
                  MYK = 0
                  MIP = CUOPS%MIP
               ELSE IF ( PSTYPE.EQ.'U' ) THEN
                  MIP = CUOPS%MIP
                  MXK = MIP
                  MYK = 1
               END IF

               IF (.NOT.ALLOCATED(XC)) ALLOCATE(XC(MIP))
               IF (.NOT.ALLOCATED(YC)) ALLOCATE(YC(MIP))
               IF (.NOT.ALLOCATED(XP)) ALLOCATE(XP(MIP))
               IF (.NOT.ALLOCATED(YP)) ALLOCATE(YP(MIP))

!              --- compute coordinates of output points                   40.51

               ITMP1  = MXC
               ITMP2  = MYC
               ITMP3  = MCGRD
               ITMP5  = MXF
               ITMP6  = MYF
               RTMP1  = XCLMIN
               RTMP2  = XCLMAX
               RTMP3  = YCLMIN
               RTMP4  = YCLMAX
               MXC    = MXCGL
               MYC    = MYCGL
               MCGRD  = MCGRDGL
               MXF    = 1
               MYF    = 1
               XCLMIN = XCGMIN
               XCLMAX = XCGMAX
               YCLMIN = YCGMIN
               YCLMAX = YCGMAX
               CALL SwashCoorOutp               (PSTYPE              ,
     &                      CUOPS%OPI           ,CUOPS%OPR           ,
     &                      CUOPS%XP            ,CUOPS%YP            ,
     &                      MIP                 ,XP                  ,
     &                      YP                  ,XC                  ,
     &                      YC                  ,KGRPGL              ,
     &                      XGRDGL              ,YGRDGL              )
               MXC    = ITMP1
               MYC    = ITMP2
               MCGRD  = ITMP3
               MXF    = ITMP5
               MYF    = ITMP6
               XCLMIN = RTMP1
               XCLMAX = RTMP2
               YCLMIN = RTMP3
               YCLMAX = RTMP4

!              --- correct problem coordinates with offset values

               DO IP = 1, MIP
                  RTMP1 = XP(IP)
                  RTMP2 = YP(IP)
                  IF (.NOT.EQREAL(RTMP1,OVEXCV(1))) XP(IP)=RTMP1+XOFFS
                  IF (.NOT.EQREAL(RTMP2,OVEXCV(2))) YP(IP)=RTMP2+YOFFS
               END DO

!              --- rewrite table output by means of collection of
!                  output locations

               IF ( RTYPE(1:3).EQ.'TAB' ) THEN
                  CALL SWCOLTAB ( RTYPE, CORQ%OQI, CORQ%IVTYP, MIP, IRQ,
     &                            BLKND, XC, YC, XP, YP )                 40.51 40.41
                  IF (STPNOW()) RETURN
               END IF

!              --- rewrite block output by means of collection of process
!                  output data

               IF ( RTYPE(1:3).EQ.'BLK' ) THEN
                  CALL SWCOLBLK ( RTYPE, CORQ%OQI, CORQ%IVTYP, CORQ%FAC,
     &                            SNAMPF, MXK, MYK, IRQ, BLKND, XC, YC )  40.51 40.41
                  IF (STPNOW()) RETURN
               END IF

               IF (ALLOCATED(XP)) DEALLOCATE(XP)                          40.51
               IF (ALLOCATED(YP)) DEALLOCATE(YP)                          40.51
               IF (ALLOCATED(XC)) DEALLOCATE(XC)                          40.51
               IF (ALLOCATED(YC)) DEALLOCATE(YC)                          40.51

  50           CONTINUE
               CORQ => CORQ%NEXTORQ

 100        CONTINUE

            IF ( NSTATC.EQ.1.AND.IT.LT.IT1 ) TIMCO = TIMCO + DT

         END DO

      END DO

!     --- close all files and delete process files

      DO IUNIT = HIOPEN+1, HIOPEN+NREOQ
        INQUIRE ( UNIT=IUNIT, OPENED=OPENED )
        IF (OPENED) CLOSE(IUNIT)
      END DO
      DO IUNIT = HIOPEN+NREOQ+1, HIOPEN+NREOQ*(NPROC+1)
        INQUIRE ( UNIT=IUNIT, OPENED=OPENED )
        IF (OPENED) CLOSE ( UNIT=IUNIT, STATUS='delete' )
      END DO

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOLTAB ( RTYPE, OQI, IVTYP, MIP, IRQ, BLKND,           40.51
     &                      XC   , YC , XP   , YP )                       40.51
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE SwashCommdata3                                                  40.41
      USE SwashCommdata4                                                  40.41
      USE SwashTimecomm
      USE outp_data
      USE M_PARALL                                                        40.31
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
!     40.51: Agnieszka Herman
!     40.51: Marcel Zijlema
!     SWASH: Marcel Zijlema
!
!  1. Updates
!
!     40.30, May  03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: further optimization
!     40.51, Feb. 05: re-design output process in parallel mode
!     SWASH, Sep. 10: extension with layer-dependent quantities
!
!  2. Purpose
!
!     Printing of table output based on point set by means of
!     collecting individual process output files
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain
!     IRQ         request number
!     IVTYP       type of variable output
!     MIP         total number of output points
!     OQI         array containing output request data
!     RTYPE       type of output request
!     XC          computational grid x-coordinate of output point
!     XP          user x-coordinate of output point
!     YC          computational grid y-coordinate of output point
!     YP          user y-coordinate of output point
!
      INTEGER   IRQ, MIP, OQI(4), IVTYP(OQI(3))
      REAL      BLKND(MXCGL,MYCGL), XC(MIP), YC(MIP), XP(MIP), YP(MIP)    40.51
      CHARACTER RTYPE*4
!
!  6. Local variables
!
!     EXIST :     logical whether a file exist or not
!     FSTR  :     an auxiliary string
!     IENT  :     number of entries
!     IF    :     first non-character in string
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPROC :     loop counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IVTYPE:     type of output quantity
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     K     :     loop counter
!     LFIELD:     actual length of a part of field OUTLIN
!     MSGSTR:     string to pass message to call MSGERR
!     NLINES:     number of lines in heading
!     NREF  :     unit reference number
!     NUMDEC:     number of decimals in the table
!     NVAR  :     number of output variables
!     OPENED:     logical whether a file is open or not
!     OUTLIN:     output line
!
      INTEGER       IENT, IF, IL, ILPOS, IP, IPROC, IUNIT,
     &              IUT, IVTYPE, IXK, IYK, JVAR, K, LFIELD, NLINES,
     &              NREF, NUMDEC, NVAR
      LOGICAL       EXIST, OPENED
      CHARACTER*80  MSGSTR
      CHARACTER*18  FSTR
      CHARACTER*512 OUTLIN
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!     TXPBLA           Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWCOLOUT
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
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           count lines of heading
!
!           open individual process output files and
!           write heading to generic output file
!
!     read output data from the proper process file and write             40.51
!     it in appropriate manner to generic output file                     40.51
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLTAB')

      NREF   = OQI(1)
      NVAR   = OQI(3)
      NLINES = 0

!     --- if generic output file exist

      IF ( NREF.NE.0 .AND. NREF.NE.PRINTF ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         FILENM = FILENM(1:ILPOS-1)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

!           --- count lines of heading

            IF ( RTYPE.NE.'TABD' ) THEN
               IF ( RTYPE.EQ.'TABP' .OR. RTYPE.EQ.'TABI' ) THEN
                  NLINES = NLINES + 7
               ELSE IF ( RTYPE.EQ.'TABT' .OR. RTYPE.EQ.'TABS' .OR.
     &                   RTYPE.EQ.'TABV' .OR. RTYPE.EQ.'TABU' ) THEN
                  NLINES = NLINES + 5
                  IF ( RTYPE.EQ.'TABT' .OR. RTYPE.EQ.'TABV' ) THEN
                     NLINES = NLINES + 1
                  ELSE
                     IF ( NSTATM.EQ.1 ) NLINES = NLINES + 2
                     NLINES = NLINES + 2 + MIP
                  END IF
                  DO JVAR = 1, NVAR
                     IVTYPE = IVTYP(JVAR)
                     IF ( OVSVTY(IVTYPE).LE.2 ) THEN
                        NLINES = NLINES + 4
                     ELSE
                        NLINES = NLINES + 8
                    END IF
                  END DO
               END IF
            END IF

!           --- open individual process output files and
!               write heading to generic output file

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            DO IPROC = 1, NPROC
               WRITE(FILENM(ILPOS-3:ILPOS),100) IPROC
 100           FORMAT('-',I3.3)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM )
                  ELSE
                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
               DO IP = 1, NLINES
                  READ (IUNIT,'(A)') OUTLIN
                  CALL TXPBLA(OUTLIN,IF,IL)
                  IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
               END DO
            END DO

         END IF

      END IF

!     --- read output data from the proper process file and write         40.51
!         it in appropriate manner to generic output file                 40.51

      IF ( NREF.NE.PRINTF ) THEN
         IF ( (RTYPE.EQ.'TABS' .OR. RTYPE.EQ.'TABU') .AND.
     &        NSTATM.EQ.1 ) THEN
            DO IPROC = 1, NPROC
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               READ (IUNIT,'(A)') OUTLIN
               CALL TXPBLA(OUTLIN,IF,IL)
               IF (IPROC.EQ.1) WRITE (NREF, '(A)') OUTLIN(1:IL)
            END DO
         END IF
         IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
            IPLOOP1 : DO IP = 1, MIP                                      40.51
               IXK = NINT(XC(IP)+100.)-99                                 41.07 40.51
               IYK = NINT(YC(IP)+100.)-99                                 41.07 40.51
               IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.         41.07 40.51
     &              IYK.GT.MYCGL ) THEN                                   41.07 40.51
                  CALL WREXCV                                             40.51
                  CYCLE IPLOOP1                                           40.51
               END IF                                                     40.51
               IPROC = 1                                                  40.51
               PROCLOOP1 : DO                                             40.51
                 IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN                40.51
                    IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC      40.51
                    READ (IUNIT,'(A)') OUTLIN                             40.51
                    CALL TXPBLA(OUTLIN,IF,IL)                             40.51
                    WRITE (NREF, '(A)') OUTLIN(1:IL)                      40.51
                    EXIT PROCLOOP1                                        40.51
                 ELSE                                                     40.51
                    IPROC = IPROC + 1                                     40.51
                    IF ( IPROC.LE.NPROC ) THEN                            40.51
                       CYCLE PROCLOOP1                                    40.51
                    ELSE                                                  40.51
                       CALL WREXCV                                        40.51
                       EXIT PROCLOOP1                                     40.51
                    END IF                                                40.51
                 END IF                                                   40.51
               END DO PROCLOOP1                                           40.51
            END DO IPLOOP1                                                40.51
         ELSE
            IPLOOP2 : DO IP = 1, MIP
               IXK = NINT(XC(IP)+100.)-99
               IYK = NINT(YC(IP)+100.)-99
               IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &              IYK.GT.MYCGL ) THEN
                  CALL WREXCV
                  CYCLE IPLOOP2
               END IF
               IPROC = 1
               PROCLOOP2 : DO
                 IF ( NINT(BLKND(IXK,IYK)).EQ.IPROC ) THEN
                    IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
                    DO JVAR = 1, NVAR
                       IVTYPE = IVTYP(JVAR)
                       IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 ) THEN
                          READ (IUNIT,'(A)') OUTLIN
                          CALL TXPBLA(OUTLIN,IF,IL)
                          WRITE (NREF, '(A)') OUTLIN(1:IL)
                          IF (OVSVTY(IVTYPE).EQ.3) THEN
                             READ (IUNIT,'(A)') OUTLIN
                             CALL TXPBLA(OUTLIN,IF,IL)
                             WRITE (NREF, '(A)') OUTLIN(1:IL)
                          ENDIF
                       ELSE IF ( IVTYPE.LE.70 ) THEN
                          DO K = 1, CEILING(REAL(KMAX+1)/10)
                             READ (IUNIT,'(A)') OUTLIN
                             CALL TXPBLA(OUTLIN,IF,IL)
                             WRITE (NREF, '(A)') OUTLIN(1:IL)
                          ENDDO
                          IF (OVSVTY(IVTYPE).EQ.3) THEN
                             DO K = 1, CEILING(REAL(KMAX+1)/10)
                                READ (IUNIT,'(A)') OUTLIN
                                CALL TXPBLA(OUTLIN,IF,IL)
                                WRITE (NREF, '(A)') OUTLIN(1:IL)
                             ENDDO
                          ENDIF
                       ELSE
                          DO K = 1, CEILING(REAL(KMAX)/10)
                             READ (IUNIT,'(A)') OUTLIN
                             CALL TXPBLA(OUTLIN,IF,IL)
                             WRITE (NREF, '(A)') OUTLIN(1:IL)
                          ENDDO
                          IF (OVSVTY(IVTYPE).EQ.3) THEN
                             DO K = 1, CEILING(REAL(KMAX)/10)
                                READ (IUNIT,'(A)') OUTLIN
                                CALL TXPBLA(OUTLIN,IF,IL)
                                WRITE (NREF, '(A)') OUTLIN(1:IL)
                             ENDDO
                          ENDIF
                       ENDIF
                    ENDDO
                    EXIT PROCLOOP2
                 ELSE
                    IPROC = IPROC + 1
                    IF ( IPROC.LE.NPROC ) THEN
                       CYCLE PROCLOOP2
                    ELSE
                       CALL WREXCV
                       EXIT PROCLOOP2
                    END IF
                 END IF
               END DO PROCLOOP2
            END DO IPLOOP2
         ENDIF
      END IF                                                              40.51
      RETURN

      CONTAINS                                                            40.51
      SUBROUTINE WREXCV                                                   40.51
      IL = 1                                                              40.51
      OUTLIN = '    '                                                     40.51
      IF (RTYPE.EQ.'TABI') THEN                                           40.51
!        --- write point sequence number as first column                  40.51
         WRITE (OUTLIN(1:8), '(I8)') IP                                   40.51
         IL = 9                                                           40.51
         OUTLIN(IL:IL) = ' '                                              40.51
      END IF                                                              40.51
      DO JVAR = 1, NVAR                                                   40.51
         IVTYPE = IVTYP(JVAR)                                             40.51
         IF (IVTYPE.EQ.40) THEN                                           40.51
            IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
!              --- for time, 18 characters are needed                     40.51
               FSTR   = '(A18)'                                           40.51
               LFIELD = 18                                                40.51
               OUTLIN(IL:IL+LFIELD-1) = CHTIME                            40.51
            ELSE
               WRITE(NREF,'(A18)') CHTIME
            ENDIF
         ELSE                                                             40.51
            IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
              IF (RTYPE.EQ.'TABD') THEN                                   40.51
                 FSTR   = FLT_TABLE                                       40.51
                 LFIELD = FLD_TABLE                                       40.51
              ELSE                                                        40.51
                 FSTR   = '(F11.X)'                                       40.51
                 LFIELD = 11                                              40.51
                 NUMDEC = MAX (0,6-NINT(LOG10(ABS(OVHEXP(IVTYPE)))))      40.51
                 IF (NUMDEC.GT.9) NUMDEC = 9                              40.51
                 WRITE (FSTR(6:6), '(I1)') NUMDEC                         40.51
              END IF                                                      40.51
!             --- write value into OUTLIN                                 40.51
              IF (IVTYPE.EQ.1) THEN                                       40.51
                 WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) XP(IP)          40.51
              ELSE IF (IVTYPE.EQ.2) THEN                                  40.51
                 WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) YP(IP)          40.51
              ELSE                                                        40.51
                 WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) OVEXCV(IVTYPE)  40.51
              END IF                                                      40.51
              IF (OVSVTY(IVTYPE).EQ.3) THEN                               40.51
                 IL = IL + LFIELD + 1                                     40.51
!                --- write second component of a vectorial quantity       40.51
                 WRITE (OUTLIN(IL:IL+LFIELD-1), FMT=FSTR) OVEXCV(IVTYPE)  40.51
              END IF                                                      40.51
            ELSE
              IF ( IVTYPE.LE.50 .OR. IVTYPE.GE.100 ) THEN
                 WRITE(NREF,FMT=FLT_SWASH) OVEXCV(IVTYPE)
                 IF (OVSVTY(IVTYPE).EQ.3) THEN
                    WRITE(NREF,FMT=FLT_SWASH) OVEXCV(IVTYPE)
                 ENDIF
              ELSE IF ( IVTYPE.LE.70 ) THEN
                 WRITE(NREF,FMT=FLT_SWASH) (OVEXCV(IVTYPE), K=0,KMAX)
                 IF (OVSVTY(IVTYPE).EQ.3) THEN
                    WRITE(NREF,FMT=FLT_SWASH) (OVEXCV(IVTYPE), K=0,KMAX)
                 ENDIF
              ELSE
                 WRITE(NREF,FMT=FLT_SWASH) (OVEXCV(IVTYPE), K=1,KMAX)
                 IF (OVSVTY(IVTYPE).EQ.3) THEN
                    WRITE(NREF,FMT=FLT_SWASH) (OVEXCV(IVTYPE), K=1,KMAX)
                 ENDIF
              ENDIF
            ENDIF
         END IF                                                           40.51
         IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
            IL = IL + LFIELD + 1                                          40.51
            OUTLIN(IL-1:IL) = '  '                                        40.51
         ENDIF
      END DO                                                              40.51
      IF (RTYPE.NE.'TABU' .AND. RTYPE.NE.'TABV') THEN
         CALL TXPBLA(OUTLIN,IF,IL)                                        40.51
         WRITE (NREF, '(A)') OUTLIN(1:IL)                                 40.51
      ENDIF
      RETURN                                                              40.51
      END SUBROUTINE WREXCV                                               40.51

      END
!****************************************************************
!
      SUBROUTINE SWCOLBLK ( RTYPE, OQI, IVTYP, FAC  , PSNAME,
     &                      MXK  , MYK, IRQ  , BLKND, XC    ,             40.51
     &                      YC   )                                        40.51
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE SwashCommdata3, ONLY: KMAX
      USE SwashTimecomm
      USE outp_data
      USE M_PARALL
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
!     40.51: Marcel Zijlema
!
!  1. Updates
!
!     40.31, Dec. 03: New subroutine
!     40.41, Jun. 04: some improvements with respect to MATLAB
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Dec. 04: optimization output with respect to COMPGRID
!     40.51, Feb. 05: further optimization
!     40.51, Feb. 05: re-design output process in parallel mode
!
!  2. Purpose
!
!     Writing of block output by means of collecting
!     individual process output files
!
!  4. Argument variables
!
!     BLKND       collected array giving node number per subdomain
!     FAC         factors of multiplication of block output
!     IRQ         request number
!     IVTYP       type of variable output
!     MXK         number of points in x-direction of output frame
!     MYK         number of points in y-direction of output frame
!     OQI         array containing output request data
!     PSNAME      name of output locations
!     RTYPE       type of output request
!     XC          computational grid x-coordinate of output point
!     YC          computational grid y-coordinate of output point
!
      INTEGER   MXK, MYK, IRQ, OQI(4), IVTYP(OQI(3))
      REAL      BLKND(MXCGL,MYCGL), XC(MXK*MYK), YC(MXK*MYK)              40.51
      REAL      FAC(OQI(3))
      CHARACTER RTYPE*4, PSNAME*8
!
!  6. Local variables
!
!     CRFORM:     format string for representing layer number
!     CTIM  :     string representing date of computation
!     DFAC  :     multiplication factor of block output
!     EXIST :     logical whether a file exist or not
!     FMAX  :     auxiliary real
!     FTIP1 :     auxiliary real
!     FTIP2 :     auxiliary real
!     IBN   :     actual node number
!     IDLA  :     lay-out indicator
!     IF    :     first non-character in string
!     IFAC  :     auxiliary integer
!     IENT  :     number of entries
!     IL    :     last non-character in string
!     ILPOS :     actual length of filename
!     IP    :     loop counter
!     IPD   :     switch for printing on paper or writing to file
!     IPROC :     loop counter
!     IREC  :     direct access file record counter
!     IUNIT :     counter for file unit numbers
!     IUT   :     auxiliary integer representing reference number
!     IVTYPE:     type of output quantity
!     IXK   :     loop counter
!     IYK   :     loop counter
!     JVAR  :     loop counter
!     K     :     counter
!     KFML  :     format length
!     KNUM  :     layer number
!     MATLAB:     indicates whether binary Matlab files are used
!     MSGSTR:     string to pass message to call MSGERR
!     NAMVAR:     name of MATLAB variable
!     NREF  :     unit reference number
!     NVAR  :     number of output variables
!     OPENED:     logical whether a file is open or not
!     VOQ   :     collected output variables
!     VOQK  :     collected layer-dependent output variables
!
      INTEGER      IBN, IDLA, IENT, IF, IFAC, IL, ILPOS, IP, IPD,
     &             IPROC, IUNIT, IUT, IXK, IYK, IVTYPE, JVAR, K, KFML,
     &             NREF, NVAR
      REAL         DFAC, FMAX, FTIP1, FTIP2
      LOGICAL      EXIST, OPENED
      INTEGER, SAVE :: IREC(MAX_OUTP_REQ)=0                               40.51 40.41
      LOGICAL, SAVE :: MATLAB=.FALSE.                                     40.41
      LOGICAL IVSTAT1, IVSTAT2, IVSTAT3
      CHARACTER*80 MSGSTR
      CHARACTER (LEN=20) :: CTIM                                          40.41
      CHARACTER (LEN=30) :: NAMVAR                                        40.41
      CHARACTER (LEN=4)  :: KNUM
      CHARACTER (LEN=6)  :: CRFORM
      REAL, ALLOCATABLE :: VOQ(:,:), VOQK(:,:,:)
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     SBLKPT           Writes block output to an ASCII file
!     STRACE           Tracing routine for debugging
!     SWRMAT           Writes block output to a binary Matlab file
!     TABHED           Prints heading
!     TXPBLA           Removes leading and trailing blanks in string      40.41
!
!  9. Subroutines calling
!
!     SWCOLOUT
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
!     if generic output file exist
!
!        if it is not open or already open thru another request
!
!           open generic output file or reset reference number
!
!           open individual process output files
!
!     read data from the proper process file and write                    40.51
!     it in appropriate manner to generic output file                     40.51
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOLBLK')

      NREF = OQI(1)
      NVAR = OQI(3)
      IDLA = OQI(4)

      IF ( RTYPE.EQ.'BLKP' ) THEN
        IPD = 1
        IF (NREF.EQ.PRINTF) CALL TABHED ('SWASH', PRINTF)
      ELSE IF ( RTYPE.EQ.'BLKD' ) THEN
        IPD = 2
      ELSE
        IPD = 3
      ENDIF

!     --- if generic output file exist

      IF ( NREF.NE.0 .AND. NREF.NE.PRINTF ) THEN

         FILENM = OUTP_FILES(OQI(2))
         ILPOS  = INDEX(FILENM, '-0')
         FILENM = FILENM(1:ILPOS-1)
         INQUIRE ( FILE=FILENM, OPENED=OPENED, NUMBER=IUT )
         MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.                       40.41
     &            INDEX (FILENM, '.mat' ).NE.0                            40.41

!        --- if it is not open or already open thru another request

         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN

!           --- open generic output file or reset reference number

            IF ( .NOT.OPENED ) THEN
               NREF = HIOPEN + IRQ
               OPEN ( UNIT=NREF, FILE=FILENM )
            ELSE
               NREF = IUT
            END IF
            OQI(1) = NREF

            IF (MATLAB .AND. .NOT.OPENED) THEN
               CLOSE(NREF)
               OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',
     &              STATUS='REPLACE',
!MatL4     &              ACCESS='DIRECT', RECL=1)
     &              ACCESS='DIRECT', RECL=4)                              41.08
               IREC(IRQ) = 1                                              40.51
            END IF

!           --- open individual process output files

            FILENM = OUTP_FILES(OQI(2))
            ILPOS  = INDEX ( FILENM, ' ' )-1
            DO IPROC = 1, NPROC
               WRITE(FILENM(ILPOS-3:ILPOS),100) IPROC
 100           FORMAT('-',I3.3)
               IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
               INQUIRE ( FILE=FILENM, EXIST=EXIST, OPENED=OPENED )
               IF ( .NOT.OPENED ) THEN
                  IF (EXIST) THEN
                     OPEN ( UNIT=IUNIT, FILE=FILENM, FORM='UNFORMATTED')
                  ELSE
                     MSGSTR= 'file '//FILENM(1:ILPOS)//' does not exist'
                     CALL MSGERR( 4, MSGSTR )
                     RETURN
                  END IF
               END IF
            END DO

         END IF

      END IF

!     --- read data from the proper process file and write                40.51
!         it in appropriate manner to generic output file                 40.51

      CTIM = CHTIME                                                       40.41
      CALL TXPBLA(CTIM,IF,IL)                                             40.41
      IF (ITMOPT.EQ.1) THEN
         CTIM(9:9)='_'                                                    40.41
      ELSE IF (ITMOPT.EQ.7) THEN
         CTIM(7:7)='_'
      ENDIF

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

      ALLOCATE(VOQ (MXK*MYK,       2))
      ALLOCATE(VOQK(MXK*MYK,0:KMAX,2))

      DO JVAR = 1, NVAR

         IVTYPE = IVTYP(JVAR)
         DFAC   = FAC(JVAR)

         VOQ  = OVEXCV(IVTYPE)
         VOQK = OVEXCV(IVTYPE)

         IF ( IVTYPE.LE.50 ) THEN
            DO IP = 1, MXK*MYK                                            40.51
               IXK = NINT(XC(IP)+100.)-99                                 41.07 40.51
               IYK = NINT(YC(IP)+100.)-99                                 41.07 40.51
               IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &              IYK.GT.MYCGL ) THEN
                  IBN   = 0
                  IPROC = NPROC+1
               ELSE
                  IBN   = NINT(BLKND(IXK,IYK))
                  IPROC = 1
               END IF
               PROCLOOP1 : DO                                             40.51
                 IF ( IBN.EQ.IPROC ) THEN                                 40.51
                    IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC      40.51
                    READ (IUNIT) VOQ(IP,1)                                40.51
                    EXIT PROCLOOP1                                        40.51
                 ELSE                                                     40.51
                    IPROC = IPROC + 1                                     40.51
                    IF ( IPROC.LE.NPROC ) THEN                            40.51
                       CYCLE PROCLOOP1                                    40.51
                    ELSE                                                  40.51
                       EXIT PROCLOOP1                                     40.51
                    END IF                                                40.51
                 END IF                                                   40.51
               END DO PROCLOOP1                                           40.51
            END DO                                                        40.51
            IF ( OVSVTY(IVTYPE).GE.3 ) THEN
               DO IP = 1, MXK*MYK                                         40.51
                  IXK = NINT(XC(IP)+100.)-99                              41.07 40.51
                  IYK = NINT(YC(IP)+100.)-99                              41.07 40.51
                  IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &                 IYK.GT.MYCGL ) THEN
                     IBN   = 0
                     IPROC = NPROC+1
                  ELSE
                     IBN   = NINT(BLKND(IXK,IYK))
                     IPROC = 1
                  END IF
                  PROCLOOP2 : DO                                          40.51
                    IF ( IBN.EQ.IPROC ) THEN                              40.51
                       IUNIT=HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC     40.51
                       READ (IUNIT) VOQ(IP,2)                             40.51
                       EXIT PROCLOOP2                                     40.51
                    ELSE                                                  40.51
                       IPROC = IPROC + 1                                  40.51
                       IF ( IPROC.LE.NPROC ) THEN                         40.51
                          CYCLE PROCLOOP2                                 40.51
                       ELSE                                               40.51
                          EXIT PROCLOOP2                                  40.51
                       END IF                                             40.51
                    END IF                                                40.51
                  END DO PROCLOOP2                                        40.51
               END DO                                                     40.51
            END IF
         ELSEIF ( IVTYPE.LE.70 ) THEN
            DO IP = 1, MXK*MYK
               IXK = NINT(XC(IP)+100.)-99
               IYK = NINT(YC(IP)+100.)-99
               IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &              IYK.GT.MYCGL ) THEN
                  IBN   = 0
                  IPROC = NPROC+1
               ELSE
                  IBN   = NINT(BLKND(IXK,IYK))
                  IPROC = 1
               END IF
               PROCLOOP3 : DO
                 IF ( IBN.EQ.IPROC ) THEN
                    IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
                    READ (IUNIT) (VOQK(IP,K,1), K=0,KMAX)
                    EXIT PROCLOOP3
                 ELSE
                    IPROC = IPROC + 1
                    IF ( IPROC.LE.NPROC ) THEN
                       CYCLE PROCLOOP3
                    ELSE
                       EXIT PROCLOOP3
                    END IF
                 END IF
               END DO PROCLOOP3
            END DO
            IF ( OVSVTY(IVTYPE).GE.3 ) THEN
               DO IP = 1, MXK*MYK
                  IXK = NINT(XC(IP)+100.)-99
                  IYK = NINT(YC(IP)+100.)-99
                  IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &                 IYK.GT.MYCGL ) THEN
                     IBN   = 0
                     IPROC = NPROC+1
                  ELSE
                     IBN   = NINT(BLKND(IXK,IYK))
                     IPROC = 1
                  END IF
                  PROCLOOP4 : DO
                    IF ( IBN.EQ.IPROC ) THEN
                       IUNIT=HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
                       READ (IUNIT) (VOQK(IP,K,2), K=0,KMAX)
                       EXIT PROCLOOP4
                    ELSE
                       IPROC = IPROC + 1
                       IF ( IPROC.LE.NPROC ) THEN
                          CYCLE PROCLOOP4
                       ELSE
                          EXIT PROCLOOP4
                       END IF
                    END IF
                  END DO PROCLOOP4
               END DO
            END IF
         ELSE
            DO IP = 1, MXK*MYK
               IXK = NINT(XC(IP)+100.)-99
               IYK = NINT(YC(IP)+100.)-99
               IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &              IYK.GT.MYCGL ) THEN
                  IBN   = 0
                  IPROC = NPROC+1
               ELSE
                  IBN   = NINT(BLKND(IXK,IYK))
                  IPROC = 1
               END IF
               PROCLOOP5 : DO
                 IF ( IBN.EQ.IPROC ) THEN
                    IUNIT = HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
                    READ (IUNIT) (VOQK(IP,K,1), K=1,KMAX)
                    EXIT PROCLOOP5
                 ELSE
                    IPROC = IPROC + 1
                    IF ( IPROC.LE.NPROC ) THEN
                       CYCLE PROCLOOP5
                    ELSE
                       EXIT PROCLOOP5
                    END IF
                 END IF
               END DO PROCLOOP5
            END DO
            IF ( OVSVTY(IVTYPE).GE.3 ) THEN
               DO IP = 1, MXK*MYK
                  IXK = NINT(XC(IP)+100.)-99
                  IYK = NINT(YC(IP)+100.)-99
                  IF ( IXK.LT.1 .OR. IYK.LT.1 .OR. IXK.GT.MXCGL .OR.
     &                 IYK.GT.MYCGL ) THEN
                     IBN   = 0
                     IPROC = NPROC+1
                  ELSE
                     IBN   = NINT(BLKND(IXK,IYK))
                     IPROC = 1
                  END IF
                  PROCLOOP6 : DO
                    IF ( IBN.EQ.IPROC ) THEN
                       IUNIT=HIOPEN+NREOQ+(NREF-HIOPEN-1)*NPROC+IPROC
                       READ (IUNIT) (VOQK(IP,K,2), K=1,KMAX)
                       EXIT PROCLOOP6
                    ELSE
                       IPROC = IPROC + 1
                       IF ( IPROC.LE.NPROC ) THEN
                          CYCLE PROCLOOP6
                       ELSE
                          EXIT PROCLOOP6
                       END IF
                    END IF
                  END DO PROCLOOP6
               END DO
            END IF
         ENDIF

         IF ( IVTYPE.LT. 4 .OR. IVTYPE.EQ. 5 .OR. IVTYPE.GT.20 .AND.
     &        IVTYPE.LT.27 .OR. IVTYPE.GT.32 .AND. IVTYPE.LT.38 ) THEN
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
     &        IVTYPE.NE.90 .AND. IVTYPE.NE.91      ) THEN
            IVSTAT3 = .TRUE.
         ELSE
            IVSTAT3 = .FALSE.
         ENDIF

         IF ( IPD.EQ.1 ) THEN
            IF ( DFAC.LE.0. ) THEN
               IF ( OVHEXP(IVTYPE).LT.0.5E10 ) THEN
                  IFAC = INT (10.+LOG10(OVHEXP(IVTYPE))) - 13
               ELSE
                  FMAX = 1000.
                  IF ( IVTYPE.LE.50 ) THEN
                     IF (OVSVTY(IVTYPE).EQ.1) THEN
                        FMAX = MAXVAL(ABS(VOQ(:,1)))
                     ELSE IF (OVSVTY(IVTYPE).EQ.3) THEN
                        FTIP1 = MAXVAL(ABS(VOQ(:,1)))
                        FTIP2 = MAXVAL(ABS(VOQ(:,2)))
                        FMAX  = MAX(FTIP1, FTIP2)
                     ENDIF
                  ELSE
                     IF (OVSVTY(IVTYPE).EQ.1) THEN
                        FMAX = MAXVAL(ABS(VOQK(:,:,1)))
                     ELSE IF (OVSVTY(IVTYPE).EQ.3) THEN
                        FTIP1 = MAXVAL(ABS(VOQK(:,:,1)))
                        FTIP2 = MAXVAL(ABS(VOQK(:,:,2)))
                        FMAX  = MAX(FTIP1, FTIP2)
                     ENDIF
                  ENDIF
                  IFAC = INT (10.+LOG10(FMAX)) - 13
               END IF
               DFAC = 10.**IFAC
            END IF
         ELSE
           IF ( DFAC.LE.0. ) DFAC = 1.
         END IF

         IF (OVSVTY(IVTYPE) .LT. 3) THEN
            IF (MATLAB) THEN
               IF ( IVSTAT1 ) THEN                                        40.94 40.41
                  NAMVAR = OVSNAM(IVTYPE)                                 40.41
                  CALL SWRMAT( MYK, MXK, NAMVAR, VOQ(1,1), NREF,
     &                         IREC(IRQ), IDLA, OVEXCV(IVTYPE) )          40.51
               ELSE IF ( IVSTAT2 ) THEN
                  DO K = 0, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      KNUM
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,1), NREF,
     &                            IREC(IRQ), IDLA, OVEXCV(IVTYPE) )
                  ENDDO
               ELSE IF ( IVSTAT3 ) THEN
                  DO K = 1, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      KNUM
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,1), NREF,
     &                            IREC(IRQ), IDLA, OVEXCV(IVTYPE) )
                  ENDDO
               ELSE IF ( IVTYPE.LE.50 ) THEN                              40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))     40.41
                  IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                  '_'//CTIM
                  CALL SWRMAT( MYK, MXK, NAMVAR, VOQ(1,1), NREF,
     &                         IREC(IRQ), IDLA, OVEXCV(IVTYPE) )          40.51
               ELSE IF ( IVTYPE.LE.70 ) THEN
                  DO K = 0, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      KNUM
                     IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                     '_'//CTIM
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,1), NREF,
     &                            IREC(IRQ), IDLA, OVEXCV(IVTYPE) )
                  ENDDO
               ELSE
                  DO K = 1, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      KNUM
                     IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                     '_'//CTIM
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,1), NREF,
     &                            IREC(IRQ), IDLA, OVEXCV(IVTYPE) )
                  ENDDO
               END IF                                                     40.41
            ELSE
               IF ( IVTYPE.LE.50 ) THEN
                  CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                        MXK, MYK, IDLA, OVLNAM(IVTYPE), VOQ(1,1))
               ELSE IF ( IVTYPE.LE.70 ) THEN
                 DO K = 0, KMAX
                    WRITE(KNUM, FMT=CRFORM) K
                    CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                        MXK, MYK, IDLA,
     &                        OVLNAM(IVTYPE)(1:KFML)//' k='//TRIM(KNUM),
     &                        VOQK(1,K,1))
                 ENDDO
               ELSE
                 DO K = 1, KMAX
                    WRITE(KNUM, FMT=CRFORM) K
                    CALL SBLKPT(IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                        MXK, MYK, IDLA,
     &                        OVLNAM(IVTYPE)(1:KFML)//' k='//TRIM(KNUM),
     &                        VOQK(1,K,1))
                 ENDDO
               ENDIF
            END IF
         ELSE
            IF (MATLAB) THEN
               IF ( IVSTAT1 ) THEN
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                     '_x'
                  CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,1), NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE) )
               ELSE IF ( IVSTAT2 ) THEN
                  DO K = 0, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_x'
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,1),
     &                            NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE))
                  ENDDO
               ELSE IF ( IVSTAT3 ) THEN
                  DO K = 1, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_x'
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,1),
     &                            NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE))
                  ENDDO
               ELSE IF ( IVTYPE.LE.50 ) THEN                              40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_x'                                           40.41
                  IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                  '_'//CTIM
                  CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,1), NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE) )  40.51
               ELSE IF ( IVTYPE.LE.70 ) THEN
                  DO K = 0, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_x'
                     IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                     '_'//CTIM
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,1),
     &                            NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE))
                  ENDDO
               ELSE
                  DO K = 1, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_x'
                     IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                     '_'//CTIM
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,1),
     &                            NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE))
                  ENDDO
               END IF                                                     40.41
               IF ( IVSTAT1 ) THEN
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                     '_y'
                  CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,2), NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE) )
               ELSE IF ( IVSTAT2 ) THEN
                  DO K = 0, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_y'
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,2),
     &                            NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE))
                  ENDDO
               ELSE IF ( IVSTAT3 ) THEN
                  DO K = 1, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_y'
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,2),
     &                            NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE))
                  ENDDO
               ELSE IF ( IVTYPE.LE.50 ) THEN                              40.41
                  NAMVAR = OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//   40.41
     &                     '_y'                                           40.41
                  IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                  '_'//CTIM
                  CALL SWRMAT( MYK, MXK, NAMVAR,
     &                 VOQ(1,2), NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE) )  40.51
               ELSE IF ( IVTYPE.LE.70 ) THEN
                  DO K = 0, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_y'
                     IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                     '_'//CTIM
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,2),
     &                            NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE))
                  ENDDO
               ELSE
                  DO K = 1, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     NAMVAR=OVSNAM(IVTYPE)(1:LEN_TRIM(OVSNAM(IVTYPE)))//
     &                      TRIM(KNUM)//'_y'
                     IF (IL.NE.1) NAMVAR = NAMVAR(1:LEN_TRIM(NAMVAR))//
     &                                     '_'//CTIM
                     CALL SWRMAT( MYK, MXK, NAMVAR, VOQK(1,K,2),
     &                            NREF, IREC(IRQ), IDLA, OVEXCV(IVTYPE))
                  ENDDO
               END IF                                                     40.41
            ELSE
               IF ( IVTYPE.LE.50 ) THEN
                  CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                         MXK, MYK, IDLA,
     &                         OVLNAM(IVTYPE)(1:52)//'  X-comp',
     &                         VOQ(1,1) )
                  CALL SBLKPT( IPD, NREF, DFAC, PSNAME, OVUNIT(IVTYPE),
     &                         MXK, MYK, IDLA,
     &                         OVLNAM(IVTYPE)(1:52)//'  Y-comp',
     &                         VOQ(1,2) )
               ELSE IF ( IVTYPE.LE.70 ) THEN
                  DO K = 0, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     CALL SBLKPT( IPD, NREF, DFAC, PSNAME,
     &                            OVUNIT(IVTYPE), MXK, MYK, IDLA,
     &                            OVLNAM(IVTYPE)(1:KFML-8)//
     &                            ' k='//TRIM(KNUM)//'  X-comp',
     &                            VOQK(1,K,1) )
                     CALL SBLKPT( IPD, NREF, DFAC, PSNAME,
     &                            OVUNIT(IVTYPE), MXK, MYK, IDLA,
     &                            OVLNAM(IVTYPE)(1:KFML-8)//
     &                            ' k='//TRIM(KNUM)//'  Y-comp',
     &                            VOQK(1,K,2) )
                  ENDDO
               ELSE
                  DO K = 1, KMAX
                     WRITE(KNUM, FMT=CRFORM) K
                     CALL SBLKPT( IPD, NREF, DFAC, PSNAME,
     &                            OVUNIT(IVTYPE), MXK, MYK, IDLA,
     &                            OVLNAM(IVTYPE)(1:KFML-8)//
     &                            ' k='//TRIM(KNUM)//'  X-comp',
     &                            VOQK(1,K,1) )
                     CALL SBLKPT( IPD, NREF, DFAC, PSNAME,
     &                            OVUNIT(IVTYPE), MXK, MYK, IDLA,
     &                            OVLNAM(IVTYPE)(1:KFML-8)//
     &                            ' k='//TRIM(KNUM)//'  Y-comp',
     &                            VOQK(1,K,2) )
                  ENDDO
               ENDIF
            END IF
         END IF

      END DO

      IF (IPD.EQ.1 .AND. NREF.EQ.PRINTF) WRITE (PRINTF, 111)

      DEALLOCATE(VOQ,VOQK)

  111 FORMAT (///)

      RETURN
      END
