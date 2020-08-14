subroutine SwashComputStruc
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
!    1.00, March 2010: New subroutine
!
!   Purpose
!
!   Performs one full simulation step with structured grid
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashComputStruc')
    !
    ! compute flow
    !
!TIMG    call SWTSTA(51)
    if ( .not.momskip ) call SwashComputFlow
!TIMG    call SWTSTO(51)
    if (STPNOW()) return
    !
    ! update water depths
    !
!TIMG    call SWTSTA(52)
    if ( kmax == 1 ) then
       call SwashUpdateDepths (   u1,   v1 )
    else
       call SwashUpdateDepths ( udep, vdep )
    endif
!TIMG    call SWTSTO(52)
    if (STPNOW()) return
    !
    ! update layer interfaces
    !
!TIMG    call SWTSTA(53)
    if ( kmax > 1 ) call SwashLayerIntfaces
!TIMG    call SWTSTO(53)
    !
    ! update mask arrays for wetting and drying
    !
!TIMG    call SWTSTA(54)
    call SwashDryWet
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! update mask arrays for free surface / pressurized flow
    !
!TIMG    call SWTSTA(54)
    if ( ifloat /= 0 ) call SwashPresFlow
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! update mask array for wave breaking
    !
!TIMG    call SWTSTA(54)
    if ( isurf /= 0 ) call SwashBreakPoint
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! compute transport
    !
!TIMG    call SWTSTA(101)
    if ( itrans /= 0 ) call SwashComputTrans
!TIMG    call SWTSTO(101)
    if (STPNOW()) return
    !
    ! compute 3D turbulence
    !
!TIMG    call SWTSTA(102)
    if ( iturb /= 0 ) call SwashComputTurb
!TIMG    call SWTSTO(102)
    if (STPNOW()) return
    !
end subroutine SwashComputStruc
