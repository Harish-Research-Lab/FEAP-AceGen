!$Id:$
      subroutine SB_activedofs(du,dofu, nelu, npde, ix, ndf)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2021: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    01/12/2020
!-----[--+---------+---------+---------+---------+---------+---------+-]
!     Purpose: Set active dof's

!     Inputs:
!         du(npde)          - dofs/pde
!         dofu(10,*)        - dofs in each pde
!         nelu(npde)        - nodes/pde
!         npde              - Number of PDE's
!         ndf               - Number of maximum FEAP dof's

!     Outputs:
!         ix(ndf,*)         - Map of active dof's
!-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'debugs.h'

!     Argument variables
      integer    :: npde, ndf
      integer    :: du(npde)
      integer    :: dofu(10,*)
      integer    :: nelu(npde)
      integer    :: ix(ndf,*)

!     Local variables
      integer    :: n,a,i

      save

!     Zero used part of array
      do n = 1,npde
        do a = 1,nelu(n)
          ix(:,a+1) = 0
        end do ! a
      end do ! n

!     Set active dofs for element
      do n = 1,npde
        do a = 1,nelu(n)
          do i = 1,du(n)
            ix(dofu(n,i),a+1) = dofu(n,i)
          end do ! i
        end do ! a
      end do ! n

!     Output
      if(debug) then
        a = maxval(nelu(1:npde))
        call iprint(ix,ndf,a+1,ndf,'IX_DOF MAP')
      endif

      end subroutine SB_activedofs
