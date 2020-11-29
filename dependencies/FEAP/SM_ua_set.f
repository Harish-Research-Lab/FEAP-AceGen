!$Id:$
      subroutine SM_ua_set(ul,ndf,nen, ua, dg, nel)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2020: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    20/11/2020
!-----[--.----+----.----+----.-----------------------------------------]
!     Purpose: Map FEAP ul(ndf,nen,*) to AceGen ua(6,*) using dg(*)

!     Inputs:
!       ul(ndf,nen,*)  - FEAP local solution values
!       ndf            - Maximum number dofs/node
!       nen            - Maxumum number nodes/element
!       dg(*)          - AceGen dofs/node

!     Outputs:
!       ua(6,nel)      - AceGen local solution parameters
!                        1: u_n+1; 2: u_n+1 - u_n; 3: du_n+1;
!                        4: v_n+1; 5: a_n+1; 6: v_n (??)
!       nel            - Number of element nodes
!-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'debugs.h'

!     Argument variables
      integer       :: ndf,nen, nel
      integer       :: dg(*)
      real (kind=8) :: ul(ndf,nen,*)
      real (kind=8) :: ua(6,*)

!     Local variables
      integer       :: a,b,i,j

!     Remap dofs from FEAP to AceGen
      b = 0
      do a = 1,nel
        do j = 1,dg(a)
          b = b + 1
          do i = 1,6
            ua(i,b) = ul(j,a,i)
          end do ! i
        end do ! j
      end do ! a

!     Output remap if debug flag set
      if(debug) then
        do i = 1,6
          call mprint(ul(1,1,i),ndf,nel,ndf,'UL')
        end do ! i
        call mprint(ua,6,b,6,'UA')
      endif

      end subroutine SM_ua_set
