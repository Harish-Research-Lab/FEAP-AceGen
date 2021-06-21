!**************************************************************
!* AceGen    7.205 MacOSX (15 Jan 21)                         *
!*           Co. J. Korelc  2020           20 Jun 21 23:00:51 *
!**************************************************************
! User     : Full professional version
! Notebook : linelasticSmall
! Evaluation time                 : 8 s     Mode  : Optimal
! Number of formulae              : 195     Method: Automatic
! Subroutine                      : elmt15_ISW01 size: 57
! Subroutine                      : elmt15_ISW03 size: 3187
! Subroutine                      : elmt15_ISW05 size: 450
! Subroutine                      : elmt15_ISW06 size: 878
! Subroutine                      : elmt15_ISW09 size: 450
! Total size of Mathematica  code : 5022 subexpressions
! Total size of Fortran code      : 17153 bytes

      subroutine elmt15 (d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2021: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!-----[--.----+----.----+----.-----------------------------------------]
!     Purpose:

!     Inputs:
!       d(*)          - Material set parameters
!       ul(ndf,nen,*) - Solution parameters
!       xl(ndm,nen)   - Element nodal coordinates
!       ix(nen1)      - Element global node numbers
!       tl(nen)       - Element vector (e.g., nodal temperatures)
!       ndf           - Maximum no dof's/node
!       ndm           - Mesh spatial dimension
!       nst           - Element matrix/vector dimension
!       isw           - Switch parameter (set by input commands)

!     Outputs:
!       s(nst,nst,2)  - Element matrix (stiffness, mass, etc.)
!       p(nst,2)      - Element vector (residual, lump mass, etc.)
!-----[--.----+----.----+----.-----------------------------------------]

      implicit none

      include 'debugs.h'
      include 'bdata.h'
      include 'cdata.h'
      include 'eldata.h'
      include 'eltran.h'
      include 'hdata.h'
      include 'iofile.h'
      include 'umac1.h'                 ! utx(1)
      include 'comblk.h'

      integer       ::  ndf             ! Max DOF's/node
      integer       ::  ndm             ! Mesh spatial dimention
      integer       ::  nst             ! Matrix 's' dimension
      integer       ::  isw             ! Switch for solution steps
      integer       ::  ix(*)           ! Element global nodes
      real (kind=8) ::  d(*)            ! Material set parameters
      real (kind=8) ::  ul(ndf,nen,*)   ! Local nodal solutions
      real (kind=8) ::  xl(ndm,*)       ! Local nodal coordinates
      real (kind=8) ::  tl(*)           ! Local nodal load array
      real (kind=8) ::  s(nst,nst,*)    ! Element matrix
      real (kind=8) ::  p(ndf,*)        ! Element vector    
      logical       ::  symmetric       ! Is matrix symmetric
      logical       ::  errck           ! For error in inputs

      character (len=50) :: datades(2)
      logical       ::  pinput
      integer       ::  dofacegen,ngpo
      real (kind=8) ::  td(10)
      real (kind=8) ::  v(1003)   ! AceGen Storage
      real (kind=8), allocatable :: ua(:,:)    ! AceGen Storage
      real (kind=8), allocatable :: k(:,:)     ! AceGen stiffness
      real (kind=8), allocatable :: c(:,:)     ! AceGen damping
      real (kind=8), allocatable :: m(:,:)     ! AceGen mass
      real (kind=8), allocatable :: r(:)       ! AceGen residual

      integer       :: i,j,kk,n                 ! Loop variables
      real (kind=8) :: gp(4,25) 

!     Limits: 10 PDE; 30 DOF total
      integer       :: npde, nelu(10),du(10),dofu(10,30)
      integer       :: hist1, hist3

      save

      symmetric=.false.

      if(isw.ne.1) then

!       Read number of PDEs
        npde = nint(d(100))

      endif
 
      if(isw.lt.0) then
        utx(1) = "elmt15"

      elseif(isw.eq.1) then                ! Input material set data
        pstyp = ndm                        ! Sets plot dimension (1,2,3)

!       Read number pdes 
        nelu(:)   = 0
        du(:)     = 0
        dofu(:,:) = 0
        call elmt15_ISW01(v,npde,du,nelu,hist1,hist3)

!       Build 'dofu'
        kk = 0
        do i = 1,npde
          do j = 1, du(i)
            kk = kk + 1
            dofu(i,j) = kk
          end do ! j
        end do ! i

!       Initialize number of history variables
        nh1 = hist1
        nh3 = hist3

!       Save number of pde
        d(100) = dble(npde)

!       Build size of AceGen arrays
        dofacegen = 0
        do n = 1,npde
          dofacegen = dofacegen + nelu(n)*du(n)
        end do ! n

        if( .not.allocated(r) ) then
          if(debug) write(*,*) ' --> Allocate: DOFACEGEN =',dofacegen
          allocate( r(dofacegen) )
          allocate( m(dofacegen,dofacegen) )
          allocate( c(dofacegen,dofacegen) )
          allocate( k(dofacegen,dofacegen) )
          allocate( ua(6,dofacegen) )
        endif

!       Material properties
        errck = pinput(td,2)
        d(1:2) = td(1:2)

!       Output material set data
        
        datades(1)="Em-Youngs modulus"

        datades(2)="[Nu]-Poisson ratio"

        write(iow,"(10x,f15.5,A3,A)")
     #     (d(i)," = ",datades(i),i=1,2)

!       Set active dof's
        call SB_activedofs(du,dofu,nelu,npde,ix,ndf)

!       Set quadrature points
        nel = maxval(nelu(1:npde))
        call SB_Int(ndm, nel, ngpo,gp)

      elseif(isw.eq.3 .or.isw.eq.6) then   ! Compute residual/tangent

!       Map dofs
        call SB_ua_set(ul,ndf,nen,ua,du,dofu,nelu,npde)

        if(isw.eq.3) then

!         Initialize element arrays
          k(:,:) = 0.0d0
          c(:,:) = 0.0d0
          m(:,:) = 0.0d0
          r(:)   = 0.0d0

!         Call the routine for ISW=3
          call elmt15_ISW03(v,d,xl,ua,k,c,m,r,hr(nh1),hr(nh2),gp,ngpo)

!         Print values for debug mode
          if(debug) then
            call mprint(k,dofacegen,dofacegen,dofacegen,'ACE_K')
            call mprint(c,dofacegen,dofacegen,dofacegen,'ACE_C')
            call mprint(m,dofacegen,dofacegen,dofacegen,'ACE_M')
            call mprint(r,1,dofacegen,1,'ACE_R')
          endif

!         Combine tangent
          k(:,:) = k(:,:)*ctan(1) + c(:,:)*ctan(2) + m(:,:)*ctan(3)

!         Fill symmetric part
          if(symmetric) then
            do j = 1,dofacegen
              do i = j+1,dofacegen
                k(i,j) = k(j,i)
              end do ! i
            end do ! j
          end if

!         Map to FEAP storage
          call SB_ace2nst(k,r,dofacegen,s,p,nst,du,dofu,nelu,npde,isw)

        elseif(isw.eq.6) then

!         Initialize element arrays
          r(:)   = 0.0d0

!         Call the routine for ISW=6
          call elmt15_ISW06(v,d,xl,ua,r,hr(nh1),hr(nh2),gp,ngpo)

!         Map to FEAP storage
          call SB_ace2nst(k,r,dofacegen,s,p,nst,du,dofu,nelu,npde,isw)

        endif

      elseif(isw.eq.4 .or.isw.eq.8) then   ! Output/plot element data

      elseif(isw.eq.5) then                ! Compute mass matrix

!       Map dofs
        call SB_ua_set(ul,ndf,nen,ua,du,dofu,nelu,npde)

!       Initialize element arrays
        m(:,:) = 0.0d0
        r(:)   = 0.0d0

!       Call the routine for ISW=5
        call elmt15_ISW05(v,d,xl,ua,m,r,hr(nh1),hr(nh2),gp,ngpo)

      elseif(isw.eq.9) then                ! Compute mass matrix

!       Map dofs
        call SB_ua_set(ul,ndf,nen,ua,du,dofu,nelu,npde)

!       Initialize element arrays
        c(:,:) = 0.0d0
        r(:)   = 0.0d0

!       Call the routine for ISW=9
        call elmt15_ISW09(v,d,xl,ua,c,r,hr(nh1),hr(nh2),gp,ngpo)

      endif

      end subroutine elmt15


!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt15_ISW01(v,npde,du,nelu,hist1,hist3)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER npde,du(10),nelu(10),hist1,hist3
      DOUBLE PRECISION v(1003)
      npde=(1)
      du(1)=(2)
      nelu(1)=(4)
      hist1=(0)
      hist3=(0)
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt15_ISW03(v,d,xl,ul,s,c,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i6
      DOUBLE PRECISION v(1003),d(2),xl(2,4),ul(6,8),s(8,8),c
     &(8,8),m(8,8),p(8),ht(0),hp(0),gp(4,4)
      v(72)=ul(1,8)
      v(71)=ul(1,7)
      v(70)=ul(1,6)
      v(69)=ul(1,5)
      v(68)=ul(1,4)
      v(67)=ul(1,3)
      v(66)=ul(1,2)
      v(65)=ul(1,1)
      v(16)=xl(2,4)
      v(15)=xl(1,4)
      v(14)=xl(2,3)
      v(875)=v(14)-v(16)
      v(13)=xl(1,3)
      v(877)=v(13)-v(15)
      v(12)=xl(2,2)
      v(871)=v(12)-v(14)
      v(11)=xl(1,2)
      v(873)=v(11)-v(13)
      v(10)=xl(2,1)
      v(874)=v(10)-v(12)
      v(870)=v(10)-v(16)
      v(9)=xl(1,1)
      v(876)=-v(11)+v(9)
      v(872)=-v(15)+v(9)
      v(8)=d(2)
      v(861)=d(1)/(1d0+v(8))
      v(163)=(v(8)*v(861))/(1d0-2d0*v(8))
      v(161)=v(861)/2d0
      v(869)=2d0*v(161)
      v(878)=v(163)+v(869)
      DO i6=1,ngpo
       v(89)=gp(1,i6)
       v(104)=((-1d0)+v(89))/4d0
       v(105)=((-1d0)-v(89))/4d0
       v(110)=v(104)*v(870)+v(105)*v(871)
       v(108)=v(104)*v(872)+v(105)*v(873)
       v(90)=gp(2,i6)
       v(106)=(1d0+v(90))/4d0
       v(103)=((-1d0)+v(90))/4d0
       v(109)=v(103)*v(874)+v(106)*v(875)
       v(107)=v(103)*v(876)+v(106)*v(877)
       v(111)=-(v(108)*v(109))+v(107)*v(110)
       v(863)=gp(4,i6)*v(111)
       v(134)=-(v(110)/v(111))
       v(148)=-(v(106)*v(134))
       v(140)=-(v(103)*v(134))
       v(135)=v(108)/v(111)
       v(151)=-(v(106)*v(135))
       v(142)=-(v(103)*v(135))
       v(136)=-(v(109)/v(111))
       v(149)=v(104)*v(136)
       v(144)=v(105)*v(136)
       v(137)=v(107)/v(111)
       v(152)=v(104)*v(137)
       v(146)=v(105)*v(137)
       v(138)=v(140)+v(149)
       v(208)=v(138)*v(161)
       v(139)=v(142)+v(152)
       v(196)=v(139)*v(163)
       v(256)=(v(138)*v(196)+v(139)*v(208))*v(863)
       v(141)=-v(140)+v(144)
       v(862)=v(141)*v(161)
       v(197)=v(141)*v(163)
       v(204)=v(197)+2d0*v(862)
       v(143)=-v(142)+v(146)
       v(864)=v(143)*v(161)
       v(198)=v(143)*v(163)
       v(271)=(v(141)*v(198)+v(143)*v(862))*v(863)
       v(258)=(v(138)*v(198)+v(139)*v(862))*v(863)
       v(216)=v(198)+2d0*v(864)
       v(265)=(v(139)*v(216)+v(138)*v(862))*v(863)
       v(264)=v(863)*(v(139)*v(197)+v(138)*v(864))
       v(257)=v(863)*(v(138)*v(204)+v(139)*v(864))
       v(145)=-v(144)+v(148)
       v(865)=v(145)*v(161)
       v(199)=v(145)*v(163)
       v(205)=v(199)+2d0*v(865)
       v(147)=-v(146)+v(151)
       v(866)=v(147)*v(161)
       v(200)=v(147)*v(163)
       v(282)=v(863)*(v(145)*v(200)+v(147)*v(865))
       v(273)=v(863)*(v(141)*v(200)+v(143)*v(865))
       v(260)=v(863)*(v(138)*v(200)+v(139)*v(865))
       v(217)=v(200)+2d0*v(866)
       v(278)=v(863)*(v(143)*v(217)+v(141)*v(865))
       v(267)=v(863)*(v(139)*v(217)+v(138)*v(865))
       v(277)=v(863)*(v(143)*v(199)+v(141)*v(866))
       v(272)=v(863)*(v(141)*v(205)+v(143)*v(866))
       v(266)=v(863)*(v(139)*v(199)+v(138)*v(866))
       v(259)=v(863)*(v(138)*v(205)+v(139)*v(866))
       v(150)=-v(148)-v(149)
       v(867)=v(150)*v(161)
       v(201)=v(150)*v(163)
       v(206)=v(201)+2d0*v(867)
       v(153)=-v(151)-v(152)
       v(868)=v(153)*v(161)
       v(202)=v(153)*v(163)
       v(289)=v(863)*(v(150)*v(202)+v(153)*v(867))
       v(284)=v(863)*(v(145)*v(202)+v(147)*v(867))
       v(275)=v(863)*(v(141)*v(202)+v(143)*v(867))
       v(262)=v(863)*(v(138)*v(202)+v(139)*v(867))
       v(218)=v(202)+2d0*v(868)
       v(287)=v(863)*(v(147)*v(218)+v(145)*v(867))
       v(280)=v(863)*(v(143)*v(218)+v(141)*v(867))
       v(269)=v(863)*(v(139)*v(218)+v(138)*v(867))
       v(286)=v(863)*(v(147)*v(201)+v(145)*v(868))
       v(283)=v(863)*(v(145)*v(206)+v(147)*v(868))
       v(279)=v(863)*(v(143)*v(201)+v(141)*v(868))
       v(274)=v(863)*(v(141)*v(206)+v(143)*v(868))
       v(268)=v(863)*(v(139)*v(201)+v(138)*v(868))
       v(261)=v(863)*(v(138)*v(206)+v(139)*v(868))
       v(158)=v(138)*v(65)+v(141)*v(67)+v(145)*v(69)+v(150)*v(71)
       v(160)=v(139)*v(66)+v(143)*v(68)+v(147)*v(70)+v(153)*v(72)
       v(165)=(v(158)+v(160))*v(163)
       v(166)=v(165)+v(158)*v(869)
       v(167)=v(161)*(v(139)*v(65)+v(138)*v(66)+v(143)*v(67)+v(141)*v
     & (68)+v(147)*v(69)+v(145)*v(70)+v(153)*v(71)+v(150)*v(72))
       v(168)=v(165)+v(160)*v(869)
       p(1)=p(1)+(v(138)*v(166)+v(139)*v(167))*v(863)
       p(2)=p(2)+(v(138)*v(167)+v(139)*v(168))*v(863)
       p(3)=p(3)+(v(141)*v(166)+v(143)*v(167))*v(863)
       p(4)=p(4)+(v(141)*v(167)+v(143)*v(168))*v(863)
       p(5)=p(5)+(v(145)*v(166)+v(147)*v(167))*v(863)
       p(6)=p(6)+(v(145)*v(167)+v(147)*v(168))*v(863)
       p(7)=p(7)+(v(150)*v(166)+v(153)*v(167))*v(863)
       p(8)=p(8)+(v(150)*v(167)+v(153)*v(168))*v(863)
       s(1,1)=s(1,1)+v(863)*((v(139)*v(139))*v(161)+(v(138)*v(138))*v
     & (878))
       s(1,2)=s(1,2)+v(256)
       s(1,3)=s(1,3)+v(257)
       s(1,4)=s(1,4)+v(258)
       s(1,5)=s(1,5)+v(259)
       s(1,6)=s(1,6)+v(260)
       s(1,7)=s(1,7)+v(261)
       s(1,8)=s(1,8)+v(262)
       s(2,1)=s(2,1)+v(256)
       s(2,2)=s(2,2)+v(863)*(v(138)*v(208)+v(139)*(v(196)+v(139)*v
     & (869)))
       s(2,3)=s(2,3)+v(264)
       s(2,4)=s(2,4)+v(265)
       s(2,5)=s(2,5)+v(266)
       s(2,6)=s(2,6)+v(267)
       s(2,7)=s(2,7)+v(268)
       s(2,8)=s(2,8)+v(269)
       s(3,1)=s(3,1)+v(257)
       s(3,2)=s(3,2)+v(264)
       s(3,3)=s(3,3)+v(863)*(v(141)*v(204)+v(143)*v(864))
       s(3,4)=s(3,4)+v(271)
       s(3,5)=s(3,5)+v(272)
       s(3,6)=s(3,6)+v(273)
       s(3,7)=s(3,7)+v(274)
       s(3,8)=s(3,8)+v(275)
       s(4,1)=s(4,1)+v(258)
       s(4,2)=s(4,2)+v(265)
       s(4,3)=s(4,3)+v(271)
       s(4,4)=s(4,4)+(v(143)*v(216)+v(141)*v(862))*v(863)
       s(4,5)=s(4,5)+v(277)
       s(4,6)=s(4,6)+v(278)
       s(4,7)=s(4,7)+v(279)
       s(4,8)=s(4,8)+v(280)
       s(5,1)=s(5,1)+v(259)
       s(5,2)=s(5,2)+v(266)
       s(5,3)=s(5,3)+v(272)
       s(5,4)=s(5,4)+v(277)
       s(5,5)=s(5,5)+v(863)*(v(145)*v(205)+v(147)*v(866))
       s(5,6)=s(5,6)+v(282)
       s(5,7)=s(5,7)+v(283)
       s(5,8)=s(5,8)+v(284)
       s(6,1)=s(6,1)+v(260)
       s(6,2)=s(6,2)+v(267)
       s(6,3)=s(6,3)+v(273)
       s(6,4)=s(6,4)+v(278)
       s(6,5)=s(6,5)+v(282)
       s(6,6)=s(6,6)+v(863)*(v(147)*v(217)+v(145)*v(865))
       s(6,7)=s(6,7)+v(286)
       s(6,8)=s(6,8)+v(287)
       s(7,1)=s(7,1)+v(261)
       s(7,2)=s(7,2)+v(268)
       s(7,3)=s(7,3)+v(274)
       s(7,4)=s(7,4)+v(279)
       s(7,5)=s(7,5)+v(283)
       s(7,6)=s(7,6)+v(286)
       s(7,7)=s(7,7)+v(863)*(v(150)*v(206)+v(153)*v(868))
       s(7,8)=s(7,8)+v(289)
       s(8,1)=s(8,1)+v(262)
       s(8,2)=s(8,2)+v(269)
       s(8,3)=s(8,3)+v(275)
       s(8,4)=s(8,4)+v(280)
       s(8,5)=s(8,5)+v(284)
       s(8,6)=s(8,6)+v(287)
       s(8,7)=s(8,7)+v(289)
       s(8,8)=s(8,8)+v(863)*(v(153)*v(218)+v(150)*v(867))
       c(1,1)=0d0+c(1,1)
       c(1,2)=0d0+c(1,2)
       c(1,3)=0d0+c(1,3)
       c(1,4)=0d0+c(1,4)
       c(1,5)=0d0+c(1,5)
       c(1,6)=0d0+c(1,6)
       c(1,7)=0d0+c(1,7)
       c(1,8)=0d0+c(1,8)
       c(2,1)=0d0+c(2,1)
       c(2,2)=0d0+c(2,2)
       c(2,3)=0d0+c(2,3)
       c(2,4)=0d0+c(2,4)
       c(2,5)=0d0+c(2,5)
       c(2,6)=0d0+c(2,6)
       c(2,7)=0d0+c(2,7)
       c(2,8)=0d0+c(2,8)
       c(3,1)=0d0+c(3,1)
       c(3,2)=0d0+c(3,2)
       c(3,3)=0d0+c(3,3)
       c(3,4)=0d0+c(3,4)
       c(3,5)=0d0+c(3,5)
       c(3,6)=0d0+c(3,6)
       c(3,7)=0d0+c(3,7)
       c(3,8)=0d0+c(3,8)
       c(4,1)=0d0+c(4,1)
       c(4,2)=0d0+c(4,2)
       c(4,3)=0d0+c(4,3)
       c(4,4)=0d0+c(4,4)
       c(4,5)=0d0+c(4,5)
       c(4,6)=0d0+c(4,6)
       c(4,7)=0d0+c(4,7)
       c(4,8)=0d0+c(4,8)
       c(5,1)=0d0+c(5,1)
       c(5,2)=0d0+c(5,2)
       c(5,3)=0d0+c(5,3)
       c(5,4)=0d0+c(5,4)
       c(5,5)=0d0+c(5,5)
       c(5,6)=0d0+c(5,6)
       c(5,7)=0d0+c(5,7)
       c(5,8)=0d0+c(5,8)
       c(6,1)=0d0+c(6,1)
       c(6,2)=0d0+c(6,2)
       c(6,3)=0d0+c(6,3)
       c(6,4)=0d0+c(6,4)
       c(6,5)=0d0+c(6,5)
       c(6,6)=0d0+c(6,6)
       c(6,7)=0d0+c(6,7)
       c(6,8)=0d0+c(6,8)
       c(7,1)=0d0+c(7,1)
       c(7,2)=0d0+c(7,2)
       c(7,3)=0d0+c(7,3)
       c(7,4)=0d0+c(7,4)
       c(7,5)=0d0+c(7,5)
       c(7,6)=0d0+c(7,6)
       c(7,7)=0d0+c(7,7)
       c(7,8)=0d0+c(7,8)
       c(8,1)=0d0+c(8,1)
       c(8,2)=0d0+c(8,2)
       c(8,3)=0d0+c(8,3)
       c(8,4)=0d0+c(8,4)
       c(8,5)=0d0+c(8,5)
       c(8,6)=0d0+c(8,6)
       c(8,7)=0d0+c(8,7)
       c(8,8)=0d0+c(8,8)
       m(1,1)=0d0+m(1,1)
       m(1,2)=0d0+m(1,2)
       m(1,3)=0d0+m(1,3)
       m(1,4)=0d0+m(1,4)
       m(1,5)=0d0+m(1,5)
       m(1,6)=0d0+m(1,6)
       m(1,7)=0d0+m(1,7)
       m(1,8)=0d0+m(1,8)
       m(2,1)=0d0+m(2,1)
       m(2,2)=0d0+m(2,2)
       m(2,3)=0d0+m(2,3)
       m(2,4)=0d0+m(2,4)
       m(2,5)=0d0+m(2,5)
       m(2,6)=0d0+m(2,6)
       m(2,7)=0d0+m(2,7)
       m(2,8)=0d0+m(2,8)
       m(3,1)=0d0+m(3,1)
       m(3,2)=0d0+m(3,2)
       m(3,3)=0d0+m(3,3)
       m(3,4)=0d0+m(3,4)
       m(3,5)=0d0+m(3,5)
       m(3,6)=0d0+m(3,6)
       m(3,7)=0d0+m(3,7)
       m(3,8)=0d0+m(3,8)
       m(4,1)=0d0+m(4,1)
       m(4,2)=0d0+m(4,2)
       m(4,3)=0d0+m(4,3)
       m(4,4)=0d0+m(4,4)
       m(4,5)=0d0+m(4,5)
       m(4,6)=0d0+m(4,6)
       m(4,7)=0d0+m(4,7)
       m(4,8)=0d0+m(4,8)
       m(5,1)=0d0+m(5,1)
       m(5,2)=0d0+m(5,2)
       m(5,3)=0d0+m(5,3)
       m(5,4)=0d0+m(5,4)
       m(5,5)=0d0+m(5,5)
       m(5,6)=0d0+m(5,6)
       m(5,7)=0d0+m(5,7)
       m(5,8)=0d0+m(5,8)
       m(6,1)=0d0+m(6,1)
       m(6,2)=0d0+m(6,2)
       m(6,3)=0d0+m(6,3)
       m(6,4)=0d0+m(6,4)
       m(6,5)=0d0+m(6,5)
       m(6,6)=0d0+m(6,6)
       m(6,7)=0d0+m(6,7)
       m(6,8)=0d0+m(6,8)
       m(7,1)=0d0+m(7,1)
       m(7,2)=0d0+m(7,2)
       m(7,3)=0d0+m(7,3)
       m(7,4)=0d0+m(7,4)
       m(7,5)=0d0+m(7,5)
       m(7,6)=0d0+m(7,6)
       m(7,7)=0d0+m(7,7)
       m(7,8)=0d0+m(7,8)
       m(8,1)=0d0+m(8,1)
       m(8,2)=0d0+m(8,2)
       m(8,3)=0d0+m(8,3)
       m(8,4)=0d0+m(8,4)
       m(8,5)=0d0+m(8,5)
       m(8,6)=0d0+m(8,6)
       m(8,7)=0d0+m(8,7)
       m(8,8)=0d0+m(8,8)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt15_ISW05(v,d,xl,ul,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i294
      DOUBLE PRECISION v(1003),d(2),xl(2,4),ul(6,8),m(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i294=1,ngpo
       m(1,1)=0d0+m(1,1)
       m(1,2)=0d0+m(1,2)
       m(1,3)=0d0+m(1,3)
       m(1,4)=0d0+m(1,4)
       m(1,5)=0d0+m(1,5)
       m(1,6)=0d0+m(1,6)
       m(1,7)=0d0+m(1,7)
       m(1,8)=0d0+m(1,8)
       m(2,1)=0d0+m(2,1)
       m(2,2)=0d0+m(2,2)
       m(2,3)=0d0+m(2,3)
       m(2,4)=0d0+m(2,4)
       m(2,5)=0d0+m(2,5)
       m(2,6)=0d0+m(2,6)
       m(2,7)=0d0+m(2,7)
       m(2,8)=0d0+m(2,8)
       m(3,1)=0d0+m(3,1)
       m(3,2)=0d0+m(3,2)
       m(3,3)=0d0+m(3,3)
       m(3,4)=0d0+m(3,4)
       m(3,5)=0d0+m(3,5)
       m(3,6)=0d0+m(3,6)
       m(3,7)=0d0+m(3,7)
       m(3,8)=0d0+m(3,8)
       m(4,1)=0d0+m(4,1)
       m(4,2)=0d0+m(4,2)
       m(4,3)=0d0+m(4,3)
       m(4,4)=0d0+m(4,4)
       m(4,5)=0d0+m(4,5)
       m(4,6)=0d0+m(4,6)
       m(4,7)=0d0+m(4,7)
       m(4,8)=0d0+m(4,8)
       m(5,1)=0d0+m(5,1)
       m(5,2)=0d0+m(5,2)
       m(5,3)=0d0+m(5,3)
       m(5,4)=0d0+m(5,4)
       m(5,5)=0d0+m(5,5)
       m(5,6)=0d0+m(5,6)
       m(5,7)=0d0+m(5,7)
       m(5,8)=0d0+m(5,8)
       m(6,1)=0d0+m(6,1)
       m(6,2)=0d0+m(6,2)
       m(6,3)=0d0+m(6,3)
       m(6,4)=0d0+m(6,4)
       m(6,5)=0d0+m(6,5)
       m(6,6)=0d0+m(6,6)
       m(6,7)=0d0+m(6,7)
       m(6,8)=0d0+m(6,8)
       m(7,1)=0d0+m(7,1)
       m(7,2)=0d0+m(7,2)
       m(7,3)=0d0+m(7,3)
       m(7,4)=0d0+m(7,4)
       m(7,5)=0d0+m(7,5)
       m(7,6)=0d0+m(7,6)
       m(7,7)=0d0+m(7,7)
       m(7,8)=0d0+m(7,8)
       m(8,1)=0d0+m(8,1)
       m(8,2)=0d0+m(8,2)
       m(8,3)=0d0+m(8,3)
       m(8,4)=0d0+m(8,4)
       m(8,5)=0d0+m(8,5)
       m(8,6)=0d0+m(8,6)
       m(8,7)=0d0+m(8,7)
       m(8,8)=0d0+m(8,8)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt15_ISW06(v,d,xl,ul,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i483
      DOUBLE PRECISION v(1003),d(2),xl(2,4),ul(6,8),p(8),ht
     &(0),hp(0),gp(4,4)
      v(549)=ul(1,8)
      v(548)=ul(1,7)
      v(547)=ul(1,6)
      v(546)=ul(1,5)
      v(545)=ul(1,4)
      v(544)=ul(1,3)
      v(543)=ul(1,2)
      v(542)=ul(1,1)
      v(493)=xl(2,4)
      v(492)=xl(1,4)
      v(491)=xl(2,3)
      v(887)=v(491)-v(493)
      v(490)=xl(1,3)
      v(889)=v(490)-v(492)
      v(489)=xl(2,2)
      v(883)=v(489)-v(491)
      v(488)=xl(1,2)
      v(885)=v(488)-v(490)
      v(487)=xl(2,1)
      v(886)=v(487)-v(489)
      v(882)=v(487)-v(493)
      v(486)=xl(1,1)
      v(888)=v(486)-v(488)
      v(884)=v(486)-v(492)
      v(485)=d(2)
      v(879)=d(1)/(1d0+v(485))
      v(640)=(v(485)*v(879))/(1d0-2d0*v(485))
      v(638)=v(879)/2d0
      v(880)=2d0*v(638)
      DO i483=1,ngpo
       v(566)=gp(1,i483)
       v(581)=((-1d0)+v(566))/4d0
       v(582)=((-1d0)-v(566))/4d0
       v(587)=v(581)*v(882)+v(582)*v(883)
       v(585)=v(581)*v(884)+v(582)*v(885)
       v(567)=gp(2,i483)
       v(583)=(1d0+v(567))/4d0
       v(580)=((-1d0)+v(567))/4d0
       v(586)=v(580)*v(886)+v(583)*v(887)
       v(584)=v(580)*v(888)+v(583)*v(889)
       v(588)=-(v(585)*v(586))+v(584)*v(587)
       v(881)=gp(4,i483)*v(588)
       v(611)=-(v(587)/v(588))
       v(625)=-(v(583)*v(611))
       v(617)=-(v(580)*v(611))
       v(612)=v(585)/v(588)
       v(628)=-(v(583)*v(612))
       v(619)=-(v(580)*v(612))
       v(613)=-(v(586)/v(588))
       v(626)=v(581)*v(613)
       v(621)=v(582)*v(613)
       v(614)=v(584)/v(588)
       v(629)=v(581)*v(614)
       v(623)=v(582)*v(614)
       v(615)=v(617)+v(626)
       v(616)=v(619)+v(629)
       v(618)=-v(617)+v(621)
       v(620)=-v(619)+v(623)
       v(622)=-v(621)+v(625)
       v(624)=-v(623)+v(628)
       v(627)=-v(625)-v(626)
       v(630)=-v(628)-v(629)
       v(635)=v(542)*v(615)+v(544)*v(618)+v(546)*v(622)+v(548)*v(627)
       v(637)=v(543)*v(616)+v(545)*v(620)+v(547)*v(624)+v(549)*v(630)
       v(642)=(v(635)+v(637))*v(640)
       v(643)=v(642)+v(635)*v(880)
       v(644)=(v(543)*v(615)+v(542)*v(616)+v(545)*v(618)+v(544)*v(620
     & )+v(547)*v(622)+v(546)*v(624)+v(549)*v(627)+v(548)*v(630))*v
     & (638)
       v(645)=v(642)+v(637)*v(880)
       p(1)=p(1)+(v(615)*v(643)+v(616)*v(644))*v(881)
       p(2)=p(2)+(v(615)*v(644)+v(616)*v(645))*v(881)
       p(3)=p(3)+(v(618)*v(643)+v(620)*v(644))*v(881)
       p(4)=p(4)+(v(618)*v(644)+v(620)*v(645))*v(881)
       p(5)=p(5)+(v(622)*v(643)+v(624)*v(644))*v(881)
       p(6)=p(6)+(v(622)*v(644)+v(624)*v(645))*v(881)
       p(7)=p(7)+(v(627)*v(643)+v(630)*v(644))*v(881)
       p(8)=p(8)+(v(627)*v(644)+v(630)*v(645))*v(881)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt15_ISW09(v,d,xl,ul,c,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i672
      DOUBLE PRECISION v(1003),d(2),xl(2,4),ul(6,8),c(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i672=1,ngpo
       c(1,1)=0d0+c(1,1)
       c(1,2)=0d0+c(1,2)
       c(1,3)=0d0+c(1,3)
       c(1,4)=0d0+c(1,4)
       c(1,5)=0d0+c(1,5)
       c(1,6)=0d0+c(1,6)
       c(1,7)=0d0+c(1,7)
       c(1,8)=0d0+c(1,8)
       c(2,1)=0d0+c(2,1)
       c(2,2)=0d0+c(2,2)
       c(2,3)=0d0+c(2,3)
       c(2,4)=0d0+c(2,4)
       c(2,5)=0d0+c(2,5)
       c(2,6)=0d0+c(2,6)
       c(2,7)=0d0+c(2,7)
       c(2,8)=0d0+c(2,8)
       c(3,1)=0d0+c(3,1)
       c(3,2)=0d0+c(3,2)
       c(3,3)=0d0+c(3,3)
       c(3,4)=0d0+c(3,4)
       c(3,5)=0d0+c(3,5)
       c(3,6)=0d0+c(3,6)
       c(3,7)=0d0+c(3,7)
       c(3,8)=0d0+c(3,8)
       c(4,1)=0d0+c(4,1)
       c(4,2)=0d0+c(4,2)
       c(4,3)=0d0+c(4,3)
       c(4,4)=0d0+c(4,4)
       c(4,5)=0d0+c(4,5)
       c(4,6)=0d0+c(4,6)
       c(4,7)=0d0+c(4,7)
       c(4,8)=0d0+c(4,8)
       c(5,1)=0d0+c(5,1)
       c(5,2)=0d0+c(5,2)
       c(5,3)=0d0+c(5,3)
       c(5,4)=0d0+c(5,4)
       c(5,5)=0d0+c(5,5)
       c(5,6)=0d0+c(5,6)
       c(5,7)=0d0+c(5,7)
       c(5,8)=0d0+c(5,8)
       c(6,1)=0d0+c(6,1)
       c(6,2)=0d0+c(6,2)
       c(6,3)=0d0+c(6,3)
       c(6,4)=0d0+c(6,4)
       c(6,5)=0d0+c(6,5)
       c(6,6)=0d0+c(6,6)
       c(6,7)=0d0+c(6,7)
       c(6,8)=0d0+c(6,8)
       c(7,1)=0d0+c(7,1)
       c(7,2)=0d0+c(7,2)
       c(7,3)=0d0+c(7,3)
       c(7,4)=0d0+c(7,4)
       c(7,5)=0d0+c(7,5)
       c(7,6)=0d0+c(7,6)
       c(7,7)=0d0+c(7,7)
       c(7,8)=0d0+c(7,8)
       c(8,1)=0d0+c(8,1)
       c(8,2)=0d0+c(8,2)
       c(8,3)=0d0+c(8,3)
       c(8,4)=0d0+c(8,4)
       c(8,5)=0d0+c(8,5)
       c(8,6)=0d0+c(8,6)
       c(8,7)=0d0+c(8,7)
       c(8,8)=0d0+c(8,8)
      ENDDO
      END
