subroutine SwashDensity
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
!    1.00, May 2013: New subroutine
!
!   Purpose
!
!   Calculates density of water relative to reference density based on temperature, salinity or sediment
!
!   Method
!
!   The Eckart's formula is used here
!
!   C. Eckart
!   Properties of water, part II. The equation of state of water and sea water at low temperatures and pressures
!   Amer. J. of Sci., vol. 256, 225-240, 1958
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashFlowdata, only: rho, rp
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter over vertical layers
    integer       :: nm       ! loop counter
    !
    real*8        :: c        ! local sediment volume concentration
    real*8        :: fac      ! a factor
    real*8        :: lambda   ! polynomial function of temperature and salinity (for sea water) in Eckart's formula
    real*8        :: p0       ! another polynomial function of temperature and salinity (for sea water) in Eckart's formula
    real*8        :: rval1    ! auxiliary real
    real*8        :: rval2    ! auxiliary real
    real*8        :: rval3    ! auxiliary real
    real*8        :: s        ! local salinity
    real*8        :: t        ! local temperature
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashDensity')
    !
    ! initialize local salinity and temperature
    !
    s = dble(salw )
    t = dble(tempw)
    !
    if ( lsal > 0 .and. ltemp > 0 ) then
       !
       ! both space varying salinity and temperature are taken into account
       !
       do k = 1, kmax
          do nm = 1, mcgrd
             !
             s = dble(rp(nm,k,lsal ))
             t = dble(rp(nm,k,ltemp))
             !
             p0        =  5890d0  + (  38d0  - 375d-3*t)*t +              3d0*s
             lambda    = 17795d-1 + (1125d-2 - 745d-4*t)*t - (38d-1 + 1d-2*t)*s
             !
             rho(nm,k) = p0/(lambda+0.698*p0)
             !
          enddo
       enddo
       !
    else if ( lsal > 0 ) then
       !
       ! space varying salinity is taken into account, temperature is constant
       !
       rval1 =  5890d0  + (  38d0  - 375d-3*t)*t
       rval2 = 17795d-1 + (1125d-2 - 745d-4*t)*t
       rval3 =    38d-1 +                 1d-2*t
       !
       do k = 1, kmax
          do nm = 1, mcgrd
             !
             s = dble(rp(nm,k,lsal))
             !
             p0        = rval1 +   3d0*s
             lambda    = rval2 - rval3*s
             !
             rho(nm,k) = p0/(lambda+0.698*p0)
             !
          enddo
       enddo
       !
    else if ( ltemp > 0 ) then
       !
       ! space varying temperature is taken into account, salinity is constant
       !
       rval1 =  5890d0  +   3d0*s
       rval2 = 17795d-1 - 38d-1*s
       !
       do k = 1, kmax
          do nm = 1, mcgrd
             !
             t = dble(rp(nm,k,ltemp))
             !
             p0        = rval1 + (  38d0           - 375d-3*t)*t
             lambda    = rval2 + (1125d-2 - 1d-2*s - 745d-4*t)*t
             !
             rho(nm,k) = p0/(lambda+0.698*p0)
             !
          enddo
       enddo
       !
    else
       !
       rho = 0.
       !
    endif
    !
    ! density based on Eckart is g/ml, hence multiply by 1000 to get g/l
    ! also set density relative to reference density
    !
    if ( lsal > 0 .or. ltemp > 0 ) rho = 1000.*rho - rhow
    !
    if ( lsed > 0 .and. lmixt ) then
       !
       ! sediment is taken into account
       !
       fac = dble(rhos-rhow)
       !
       do k = 1, kmax
          do nm = 1, mcgrd
             !
             c = dble(rp(nm,k,lsed))
             !
             rho(nm,k) = c * fac + rho(nm,k) * (1d0 - c)
             !
          enddo
       enddo
       !
    endif
    !
end subroutine SwashDensity
