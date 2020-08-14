subroutine SwashAverOutp ( iquan )
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
!    1.00, March 2012: New subroutine
!
!   Purpose
!
!   Computes setup, wave height, mean current, mean constituents and mean turbulence quantities for output purposes
!
!   Method
!
!   The time-averaging procedure is employed where the requested quantities are averaged over the
!   duration in the final stage of the simulation. The simulation period is supposed to be long
!   enough to establish steady-state conditions.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
    use outp_data
    use SwashTimecomm
!
    implicit none
!
!   Argument variables
!
    integer, intent(in) :: iquan ! indicate type of quantity needed for update
                                 ! =1; surface elevation or velocity
                                 ! =2; transport constituent
                                 ! =3; turbulence quantities
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: l        ! loop counter over constituents
    integer       :: nm       ! loop counter over grid points
    !
    real*8        :: dif      ! interval between end time and current time
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashAverOutp')
    !
    dif = tfinc - timco
    !
    if ( lwavoutp .and. iquan == 1 ) then
       !
       ! setup and wave height
       !
       if ( .not. dif > twavoutp ) then
          !
          nwavoutp = nwavoutp + 1.
          !
          setup  = setup  + s1
          etavar = etavar + s1*s1
          !
       endif
       !
       if ( abs(dif) < 0.99*dt .and. nwavoutp > 1. ) then
          !
          setup = setup / nwavoutp
          !
          hsig = 4. * sqrt( max( etavar / (nwavoutp-1.) - setup*setup, 0. ) )
          !
          do nm = 1, mcgrd
             if ( .not. hsig(nm) > 0. ) setup(nm) = 0.
          enddo
          !
       endif
       !
    endif
    !
    if ( lcuroutp .and. iquan == 1 ) then
       !
       ! mean current
       !
       if ( .not. dif > tcuroutp ) then
          !
          ncuroutp = ncuroutp + 1.
          !
          mvelu = mvelu + u1
          if ( .not.oned ) mvelv = mvelv + v1
          !
       endif
       !
       if ( abs(dif) < 0.99*dt .and. ncuroutp > 0. ) then
          !
          mvelu = mvelu / ncuroutp
          if ( .not.oned ) mvelv = mvelv / ncuroutp
          !
       endif
       !
    endif
    !
    if ( ltraoutp .and. iquan == 2 ) then
       !
       ! mean constituents
       !
       if ( .not. dif > ttraoutp ) then
          !
          ntraoutp = ntraoutp + 1.
          !
          do l = 1, ltrans
             mcons(:,:,l) = mcons(:,:,l) + rp(:,:,l)
          enddo
          !
       endif
       !
       if ( abs(dif) < 0.99*dt .and. ntraoutp > 0. ) then
          !
          do l = 1, ltrans
             mcons(:,:,l) = mcons(:,:,l) / ntraoutp
          enddo
          !
       endif
       !
    endif
    !
    if ( lturoutp .and. iquan == 3 ) then
       !
       ! mean turbulence quantities
       !
       if ( .not. dif > tturoutp ) then
          !
          nturoutp = nturoutp + 1.
          !
          mtke(:,:) = mtke(:,:) + rtur(:,:,1)
          meps(:,:) = meps(:,:) + rtur(:,:,2)
          !
       endif
       !
       if ( abs(dif) < 0.99*dt .and. nturoutp > 0. ) then
          !
          mtke = mtke / nturoutp
          meps = meps / nturoutp
          !
       endif
       !
    endif
    !
end subroutine SwashAverOutp
