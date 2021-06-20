      subroutine SB_ace2nst(sall,pall,ngdof, s,p,nste,
     &                      du,dofu,nelu,npde,isw)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2021: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/12/2020
!-----[--+---------+---------+---------+---------+---------+---------+-]
!     Purpose: Move residual and tangent from AceGen to Feap storage.
!              Allows for general setting of 'ndf' to mix elements.
!              Blocked PDE form

!     Inputs:
!         sall(ngdof,ngdof) - AceGen tangent/mass array
!         pall(ngdof)       - AceGen residual/diagonal mass array
!                             N.B. Sign reversed w/r FEAP
!         ngdof             - Size of AceGen arrays
!         nste              - Size of FEAP arrays
!         npde              - Numbe pdes
!         du(npde)          - dofs/pde
!         dofu(10,*)        - dofs in each pde
!         nelu(npde)        - nodes/pde

!     Note:    sa(nel) - FEAP array pointer
!              la      - FEAP pointer to element equations NOT CODED

!     Outputs:
!         s(nste,nste)      - FEAP tangent/mass array
!         p(nste)           - FEAP residual/diagonal mass array
!-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'eldata.h'  ! nel
      include   'qudshp.h'  ! sa(*)

!     Argument variables
      integer       :: ngdof,nste, npde, isw
      integer       :: du(npde), dofu(10,*), nelu(npde)
      real (kind=8) :: sall(ngdof,ngdof), pall(ngdof)
      real (kind=8) :: s(nste,nste), p(nste)

!     Local variables
      integer       :: ar,ac, br,bc, jr,jc, nr,nc

      save

!     Repack stiffness & residual
      br = 0
      do nr = 1,npde
        do ar = 1,nelu(nr)
          do jr = 1,du(nr)
            br = br + 1
!           Assign residual
            p(sa(ar)+dofu(nr,jr)) = -pall(br)
            bc = 0
            do nc = 1,npde
              do ac = 1,nelu(nc)
                do jc = 1,du(nc)
                  bc = bc + 1
!                 Assign tangent
                  s(sa(ar)+dofu(nr,jr),sa(ac)+dofu(nc,jc)) = sall(br,bc)
                end do ! jc
              end do ! ac
            end do ! nc
          end do ! jr
        end do ! ar
      end do ! nr

      end subroutine SB_ace2nst
