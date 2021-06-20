!$Id:$
      subroutine SB_Int(ndm,nen, ngpo,gp)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2021: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    13/01/2021
!-----[--+---------+---------+---------+---------+---------+---------+-]
!     Purpose: Integration cover for AceGen - FEAP use

!     Inputs:
!       ndm             - Spatial dimension of element
!       nen             - Maximum number nodes/element

!     Outputs:
!       ngpo            - Number of quadrature points
!       gp(4,ngpo)      - Quadrature points and weights
!-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      integer       :: ndm, nen, ngpo
      real (kind=8) :: gp(4,*)

!     2-dimensional elements
      if(ndm.eq.2) then
        if(nen.eq.3) then                     ! Order 1 triangle
          call Int2dTri(3,ngpo,gp)
        elseif(nen.eq.4) then                 ! Order 1 quadrilateral
          call Int2dQuad(2,ngpo,gp)
        elseif(nen.eq.6) then                 ! Order 2 triangle
          call Int2dTri(6,ngpo,gp)
        elseif(nen.eq.10) then                ! Order 3 triangle
          call Int2dTri(7,ngpo,gp)
        elseif(nen.eq.8 .or. nen.eq.9) then   ! Order 2 quadrilateral
          call Int2dQuad(3,ngpo,gp)
        elseif(nen.eq.12 .or. nen.eq.16) then ! Order 3 quadrilateral
          call Int2dQuad(4,ngpo,gp)
        endif
      elseif(ndm.eq.3) then
        if(nen.eq.4) then                     ! Order 1 tetrahedron
          call Int3dTet(2,ngpo,gp)
        elseif(nen.eq.8) then                 ! Order 1 hexahedron
          call Int3dHex(2,ngpo,gp)
        elseif(nen.eq.10) then                ! Order 2 tetrahedron
          call Int3dTet(4,ngpo,gp)
        elseif(nen.eq.15) then                ! Order 3 tetrahedron
          call Int3dTet(6,ngpo,gp)
        elseif(nen.eq.20 .or. nen.eq.27) then ! Order 2 hexahedron
          call Int3dHex(3,ngpo,gp)
        elseif(nen.eq.64) then                ! Order 3 hexahedron
          call Int3dHex(4,ngpo,gp)
        endif
      endif

      end subroutine SB_Int
