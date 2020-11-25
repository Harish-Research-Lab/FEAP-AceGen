!$Id:$
      subroutine SM_activedofs(dg, nel, ix, ndf)
      implicit   none

!     Argument variables
      integer    :: nel, ndf
      integer    :: dg(nel)
      integer    :: ix(ndf,*)

!     Local variables
      integer    :: a,i

!     Set active dofs for element
      do a = 1,nel
        ix(:,a+1) = 0
        do i = 1,dg(a)
          ix(i,a+1) = i
        end do ! i
      end do ! a

      end subroutine SM_activedofs
