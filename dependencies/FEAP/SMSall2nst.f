      subroutine SMSall2nst(sall,pall, s,p, ngdof, nste, dg, ndg)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2020: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/11/2020
!-----[--+---------+---------+---------+---------+---------+---------+-]
!     Purpose: Move residual and tangent from AceGen to Feap storage.
!              Allows for general setting of 'ndf' to mix elements.

!     Inputs:
!         sall(ngdof,ngdof) - AceGen tangent/mass array
!         pall(ngdof)       - AceGen residual/diagonal mass array
!                             N.B. Sign reversed w/r FEAP
!         ngdof             - Size of AceGen arrays
!         nste              - Size of FEAP arrays
!         dg(ndg)           - Number DOF for each node/elmt eq.
!         ndg               - Size of AceGen dg array

!     Note:    sa(nel) - FEAP array pointer
!              la      - FEAP pointer to element equations

!     Outputs:
!         s(nste,nste)      - FEAP tangent/mass array
!         p(nste)           - FEAP residual/diagonal mass array
!-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'eldata.h'  ! nel
      include   'qudshp.h'  ! sa(*), la, ga

!     Argument variables
      integer       :: ngdof,nste, ndg
      integer       :: dg(*)
      real (kind=8) :: sall(ngdof,ngdof), pall(ngdof)
      real (kind=8) :: s(nste,nste), p(nste)

!     Local variables
      integer       :: i,j, ii,jj, i1,j1

!     Repack stiffness and residual
      j1 = 0
      do j = 1,nel
!       Residual
        do jj = 1,dg(j)
          p(sa(j)+jj) = -pall(j1+jj)
        end do ! jj
!       Tangent
        i1 = 0
        do i = 1,nel
          do jj = 1,dg(j)
            do ii = 1,dg(i)
              s(sa(i)+ii,sa(j)+jj) = sall(i1+ii,j1+jj)
            end do ! ii
          end do ! i
          i1 = i1 + dg(i)
        end do ! jj
        j1 = j1 + dg(j)
      end do ! j

!     Set element variables
      do j = nel+1,ndg
        do jj = 1,dg(j)
          p(la+jj) = -pall(j1+jj)
          i1 = 0
          do i = 1,nel
            do ii = 1,dg(i)
              s(sa(i)+ii,la+jj) = sall(i1+ii,j1+jj)
              s(la+jj,sa(i)+ii) = sall(j1+jj,i1+ii)
            end do ! ii
            i1 = i1 + dg(i)
          end do ! i
          do i = nel+1,ndg
            do ii = 1,dg(i)
              s(la+ii,la+jj) = sall(i1+ii,j1+jj)
            end do ! ii
            i1 = i1 + dg(i)
          end do ! i
        end do ! jj
        j1 = j1 + dg(j)
      end do ! j

      end subroutine SMSall2nst
