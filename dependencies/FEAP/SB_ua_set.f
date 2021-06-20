!$Id:$
      subroutine SB_ua_set(ul,ndf,nen, ua, du, dofu, nelu, npde)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2021: Regents of the University of California
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
!       du(npde,)      - AceGen dofs/node
!       dofu(10,*)     - AceGen dofs/node
!       nelu(*)        - Number of nodes/pde
!       npde           - Number of pde

!     Outputs:
!       ua(6,nel)      - AceGen local solution parameters
!                        1: u_n+1; 2: u_n+1 - u_n; 3: du_n+1;
!                        4: v_n+1; 5: a_n+1; 6: v_n (??)
!-----[--.----+----.----+----.-----------------------------------------]
      implicit   none

      include   'debugs.h'

!     Argument variables
      integer       :: ndf,nen, npde
      integer       :: du(*),dofu(10,*),nelu(*)
      real (kind=8) :: ul(ndf,nen,*)
      real (kind=8) :: ua(6,*)

!     Local variables
      integer       :: a,b,i,j, n, nel

      save

!     Remap dofs from FEAP to AceGen
      b = 0
      do n = 1,npde
        do a = 1,nelu(n)
          do j = 1,du(n)
            b = b + 1
            do i = 1,6
              ua(i,b) = ul(dofu(n,j),a,i)
            end do ! i
          end do ! j
        end do ! a
      end do ! n

!     Output remap if debug flag set
      if(debug) then
        nel = maxval(nelu(1:npde))
        do i = 1,6
          call mprint(ul(1,1,i),ndf,nel,ndf,'UL_FEAP')
        end do ! i
        call mprint(ua,6,b,6,'UA_AceGen')
      endif

      end subroutine SB_ua_set
