subroutine SwashVertVisc
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
!    1.00, February 2011: New subroutine
!
!   Purpose
!
!   Calculates vertical eddy viscosity coefficient
!
!   Method
!
!   The standard k-eps model is utilized to account for vertical mixing
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata
!
    implicit none
!
!   Parameter variables
!
    real, parameter :: cmu = 0.09    ! closure constant for standard k-eps model
    real, parameter :: eps = 1.e-8   ! small number used to prevent division by zero
    real, parameter :: maxvisc = 10. ! maximal vertical viscosity
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashVertVisc')
    !
    vnu3d(:,:) = cmu * rtur(:,:,1) * rtur(:,:,1) / ( rtur(:,:,2) + eps )
    !
    if ( iturb == 1 ) then
       !
       ! to account for other forms of unresolved mixing, background viscosity is included, if appropriate
       !
       vnu3d = max( vnu3d, bvisc )
       !
       ! add kinematic viscosity
       !
       vnu3d = vnu3d + kinvis
       !
       ! if vertical viscosity is too large, apply clipping
       !
       vnu3d = min( vnu3d, maxvisc )
       !
    endif
    !
end subroutine SwashVertVisc
