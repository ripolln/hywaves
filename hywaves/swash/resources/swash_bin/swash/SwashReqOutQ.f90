subroutine SwashReqOutQ ( found )
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
!   Authors
!
!    1.00: Marcel Zijlema
!
!   Updates
!
!    1.00, February 2010: New subroutine
!
!   Purpose
!
!   Reading and processing of the output quantities
!
!   Modules used
!
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3
    use SwashCommdata4
    use outp_data
    use m_parall
!
    implicit none
!
!   Argument variables
!
    logical, intent(inout)      :: found          ! keyword found
!
!   Local variables
!
    integer                     :: i              ! loop counter
    integer                     :: idlao          ! lay-out indicator
    integer, save               :: ient = 0       ! number of entries in this subroutine
    integer                     :: ierr           ! error indicator: ierr=0: no error, otherwise error
    integer                     :: ilpos          ! actual length of error message filename
    integer                     :: ivtype         ! type of output quantity
    integer                     :: mip            ! number of output points
    integer                     :: nref           ! reference number
    integer                     :: nvar           ! actual number of output quantities
    !
    real                        :: dfac           ! multiplication factor of block output
    !
    character(16)               :: psname         ! name of point set
    character(1)                :: pstype         ! type of point set ('F', 'H', 'U', etc.)
    character(4)                :: rtype          ! type of output request ('BLKD', 'BLKP', 'TABP', 'TABS', etc.)
    !
    logical                     :: KEYWIS         ! indicates whether keyword in user manual is found or not
    logical, save               :: lorq = .false. ! indicates whether there is request output or not
    !
    type(orqdat), pointer       :: orqtmp         ! list containing parameters for request output
    type(orqdat), save, pointer :: corq           ! current item in list of request output
    !
    type auxpt                                    ! auxiliary linking list
       integer              :: i
       real                 :: r
       type(auxpt), pointer :: nexti
    end type auxpt
    type(auxpt), target     :: frst
    type(auxpt), pointer    :: curr, tmp
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashReqOutQ')
    !
    ! ==========================================================================
    !
    ! BLOCK   'sname'  HEAD / NOHEAD  'fname' (LAY-OUT [idla])                 &
    !       < TSEC|XP|YP|DEP|BOTL|WATL|DRAF|VMAG|VDIR|VEL|VKSI|VETA|           &
    !         PRESS|NHPRES|QMAG|QDIR|DISCH|QKSI|QETA|VORT|WMAG|WDIR|WIND|      &
    !         FRC|USTAR|UFRIC|ZK|HK|VMAGK|VDIRK|VELK|VKSIK|VETAK|              &
    !         VZ|VOMEGA|QMAGK|QDIRK|DISCHK|QKSIK|QETAK|PRESSK|NHPRSK|          &
    !         TKE|EPS|VISC|HS|HRMS|SETUP|MVMAG|MVDIR|MVEL|MVKSI|MVETA|         &
    !         MVMAGK|MVDIRK|MVELK|MVKSIK|MVETAK|MTKE|MEPS|MVISC|               &
    !         SAL|TEMP|SED|MSAL|MTEMP|MSED|SALK|TEMPK|SEDK|MSALK|MTEMPK|       &
    !         MSEDK [unit] >                                                   &
    !         (OUTPUT [tbegblk] [deltblk] SEC/MIN/HR/DAY)
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('BLO') ) then
       !
       call SWNMPS ( psname, pstype, mip, ierr )
       if ( ierr /= 0 ) goto 800
       if ( pstype /= 'F' .and. pstype /= 'H' .and. pstype /= 'U' ) then
          call msgerr (2, 'set of output locations is not correct type')
          goto 800
       endif
       if ( psname(1:6) == 'NOGRID' ) then
          call msgerr (2, 'empty set is not allowed')
          goto 800
       endif
       !
       ! output frame exists
       !
       allocate(orqtmp)
       nreoq = nreoq + 1
       if ( nreoq > max_outp_req ) call msgerr (2, 'too many output requests')
       !
       idlao = 1
       !
       call INKEYW ('REQ',' ')
       if ( KEYWIS('NOHEAD') .or. KEYWIS ('FIL') ) then
          call INCSTR ('FNAME', FILENM, 'REQ', ' ')
          dfac  = 1.
          rtype = 'BLKD'
       else
          call IGNORE ('HEAD')
          call INCSTR ('FNAME', FILENM, 'STA', ' ')
          dfac  = -1.
          rtype = 'BLKP'
          if ( index (FILENM, '.MAT' ) /= 0 .or. index (FILENM, '.mat' ) /= 0 ) then
             call msgerr (4, 'no header allowed for Matlab files')
             return
          endif
       endif
       if ( FILENM == ' ' ) then
          nref = PRINTF
       else
          nref = 0
          !
          ! append node number to FILENM in case of parallel computing
          !
          if ( parll ) then
             ilpos = index ( FILENM, ' ' )-1
             write (FILENM(ilpos+1:ilpos+4),201) inode
          endif
       endif
       call INKEYW ('STA', ' ')
       if ( KEYWIS('LAY') ) then
          call ININTG ('IDLA', idlao, 'REQ', 0)
          call INKEYW ('REQ', ' ')
          if ( idlao /= 1 .and. idlao /= 3 .and. idlao /= 4 ) call msgerr (2, 'illegal value for [idla]')
       endif
       orqtmp%oqr(1) = -1.
       orqtmp%oqr(2) = -1.
       orqtmp%rqtype = rtype
       orqtmp%psname = psname
       orqtmp%oqi(1) = nref
       orqtmp%oqi(2) = nreoq
       nvar = 0
       orqtmp%oqi(3) = nvar
       orqtmp%oqi(4) = idlao
       outp_files(nreoq) = FILENM
       !
       ! read types of variables to be printed in block
       !
       frst%i = 0
       frst%r = 0.
       nullify(frst%nexti)
       curr => frst
 100   call svartp (ivtype)
       if ( ivtype == 98 ) goto 110
       if ( ivtype /= 999 ) then
          call INREAL ('UNIT', dfac, 'STA', -1.)
          if ( ovsvty(ivtype) == 7 ) then
             call msgerr (2, 'type of output not allowed for this quantity')
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          else if ( ivtype > 0 ) then
             nvar = nvar + 1
             allocate(tmp)
             tmp%i = ivtype
             tmp%r = dfac
             nullify(tmp%nexti)
             curr%nexti => tmp
             curr => tmp
             if ( ivtype == 22 .or. ivtype == 23 .or. ivtype == 24 ) lwavoutp = .true.
             if ( ivtype >= 33 .and. ivtype <= 37 ) lcuroutp = .true.
             if ( ivtype == 42 .or. ivtype == 43 .or. ivtype == 44 ) ltraoutp = .true.
             if ( ivtype == 57 .or. ivtype == 58 .or. ivtype == 59 ) lturoutp = .true.
             if ( ivtype >= 84 .and. ivtype <= 88 ) lcuroutp = .true.
             if ( ivtype == 92 .or. ivtype == 93 .or. ivtype == 94 ) ltraoutp = .true.
          endif
          goto 100
       endif
       !
 110   if ( nvar > 0 ) then
          allocate(orqtmp%ivtyp(nvar))
          allocate(orqtmp%fac(nvar))
          curr => frst%nexti
          do i = 1, nvar
             orqtmp%ivtyp(i) = curr%i
             orqtmp%fac  (i) = curr%r
             curr => curr%nexti
          enddo
          deallocate(tmp)
       endif
       !
       if ( ivtype == 98 ) then
          if ( nstatm == 0 ) call msgerr (3, 'time information not allowed in stationary mode')
          nstatm = 1
          call INCTIM (ITMOPT, 'TBEG', orqtmp%oqr(1), 'REQ', 0d0)
          call INITVD ('DELT', orqtmp%oqr(2), 'REQ', 0d0)
       endif
       !
       orqtmp%oqi(3) = nvar
       nullify(orqtmp%nextorq)
       if ( .not.lorq ) then
           forq = orqtmp
           corq => forq
           lorq = .true.
       else
          corq%nextorq => orqtmp
          corq => orqtmp
       endif
       goto 800
       !
    endif
    !
    ! ==========================================================================
    !
    ! TABLE   'sname'  HEAD | NOHEAD | STAB | SWASH | IND  'fname'             &
    !       < TIME|TSEC|XP|YP|DIST|DEP|BOTL|WATL|DRAF|VMAG|VDIR|VEL|VKSI|VETA| &
    !         PRESS|NHPRES|QMAG|QDIR|DISCH|QKSI|QETA|VORT|WMAG|WDIR|WIND|FRC|  &
    !         USTAR|UFRIC|ZK|HK|VMAGK|VDIRK|VELK|VKSIK|VETAK|VZ|VOMEGA|        &
    !         QMAGK|QDIRK|DISCHK|QKSIK|QETAK|PRESSK|NHPRSK|TKE|EPS|VISC|       &
    !         HS|HRMS|SETUP|MVMAG|MVDIR|MVEL|MVKSI|MVETA|MVMAGK|MVDIRK|        &
    !         MVELK|MVKSIK|MVETAK|MTKE|MEPS|MVISC|SAL|TEMP|SED|MSAL|MTEMP|     &
    !         MSED|SALK|TEMPK|SEDK|MSALK|MTEMPK|MSEDK|FORCEX|FORCEY|FORCEZ|    &
    !         MOMX|MOMY|MOMZ >                                                 &
    !         (OUTPUT [tbegtbl] [delttbl] SEC/MIN/HR/DAY)
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('TAB') ) then
       !
       call SWNMPS ( psname, pstype, mip, ierr )
       if ( ierr /= 0 ) goto 800
       !
       ! output points exist
       !
       allocate(orqtmp)
       nreoq = nreoq + 1
       if ( nreoq > max_outp_req ) call msgerr (2, 'too many output requests')
       !
       call INKEYW ('STA','HEAD')
       if ( KEYWIS('NOHEAD') .or. KEYWIS ('FIL') ) then
          rtype = 'TABD'
       else if ( KEYWIS ('IND') ) then
          rtype = 'TABI'
       else if ( KEYWIS ('SWASH') ) then
          rtype = 'TABS'
       else if ( KEYWIS ('STAB') ) then
          rtype = 'TABT'
       else
          call IGNORE ('HEAD')
          rtype = 'TABP'
       endif
       orqtmp%oqr(1) = -1.
       orqtmp%oqr(2) = -1.
       orqtmp%rqtype = rtype
       call INCSTR ('FNAME', FILENM, 'STA', ' ')
       if ( FILENM /= '    ' ) then
          nref = 0
          !
          ! append node number to FILENM in case of parallel computing
          !
          if ( parll ) then
             ilpos = index ( FILENM, ' ' )-1
             write (FILENM(ilpos+1:ilpos+4),201) inode
          endif
       else
          nref = PRINTF
       endif
       orqtmp%psname = psname
       orqtmp%oqi(1) = nref
       orqtmp%oqi(2) = nreoq
       outp_files(nreoq) = FILENM
       !
       nvar = 0
       orqtmp%oqi(3) = nvar
       !
       ! read types of variables to be printed in table
       !
       frst%i = 0
       nullify(frst%nexti)
       curr => frst
 120   call svartp (ivtype)
       if ( ivtype == 98 ) goto 130
       if ( ivtype /= 999 ) then
          if ( ivtype < 100 .and. ivtype /= 41 .and. psname(1:6) == 'NOGRID' ) then
             call msgerr (2, 'empty set not allowed for this quantity')
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
          endif
          if ( ivtype == 113 .and. .not.oned ) then
             call msgerr (1, 'wave runup height will not be computed in 2D mode')
          endif
          if ( ovsvty(ivtype) == 7 ) then
             call msgerr (2, 'type of output not allowed for this quantity')
             write (PRINTF, *) ' -> ', ovsnam(ivtype)
           else if ( ivtype > 0 ) then
             nvar = nvar + 1
             allocate(tmp)
             tmp%i = ivtype
             nullify(tmp%nexti)
             curr%nexti => tmp
             curr => tmp
             if ( ivtype == 22 .or. ivtype == 23 .or. ivtype == 24 ) lwavoutp = .true.
             if ( ivtype >= 33 .and. ivtype <= 37 ) lcuroutp = .true.
             if ( ivtype == 42 .or. ivtype == 43 .or. ivtype == 44 ) ltraoutp = .true.
             if ( ivtype == 57 .or. ivtype == 58 .or. ivtype == 59 ) lturoutp = .true.
             if ( ivtype >= 84 .and. ivtype <= 88 ) lcuroutp = .true.
             if ( ivtype == 92 .or. ivtype == 93 .or. ivtype == 94 ) ltraoutp = .true.
             call INKEYW ('STA', ' ')
             if ( KEYWIS('UNIT') ) then
                call msgerr (1, '[unit] is ignored in this version')
             endif
          endif
          goto 120
       endif
       !
 130   if ( nvar > 0 ) then
          allocate(orqtmp%ivtyp(nvar))
          curr => frst%nexti
          do i = 1, nvar
             orqtmp%ivtyp(i) = curr%i
             curr => curr%nexti
          enddo
          deallocate(tmp)
       endif
       allocate(orqtmp%fac(0))
       !
       if ( ivtype == 98 ) then
          if ( nstatm == 0 ) call msgerr (3, 'time information not allowed in stationary mode')
          nstatm = 1
          call INCTIM (ITMOPT, 'TBEG', orqtmp%oqr(1), 'REQ', 0d0)
          call INITVD ('DELT', orqtmp%oqr(2), 'REQ', 0d0)
       endif
       orqtmp%oqi(3) = nvar
       nullify(orqtmp%nextorq)
       if ( .not.lorq ) then
          forq = orqtmp
          corq => forq
          lorq = .true.
       else
          corq%nextorq => orqtmp
          corq => orqtmp
       endif
       goto 800
    endif
    !
    ! command not found
    !
    return
    !
 800 found = .true.
    return
    !
 201 format('-',i3.3)
    !
end subroutine SwashReqOutQ
