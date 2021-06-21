!**************************************************************
!* AceGen    7.205 MacOSX (15 Jan 21)                         *
!*           Co. J. Korelc  2020           21 Jun 21 00:04:09 *
!**************************************************************
! User     : Full professional version
! Notebook : linelasticTransient
! Evaluation time                 : 9 s     Mode  : Optimal
! Number of formulae              : 265     Method: Automatic
! Subroutine                      : elmt16_ISW01 size: 57
! Subroutine                      : elmt16_ISW03 size: 3530
! Subroutine                      : elmt16_ISW05 size: 893
! Subroutine                      : elmt16_ISW06 size: 1065
! Subroutine                      : elmt16_ISW09 size: 450
! Total size of Mathematica  code : 5995 subexpressions
! Total size of Fortran code      : 19903 bytes

      subroutine elmt16 (d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

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

      character (len=50) :: datades(3)
      logical       ::  pinput
      integer       ::  dofacegen,ngpo
      real (kind=8) ::  td(10)
      real (kind=8) ::  v(1125)   ! AceGen Storage
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
        utx(1) = "elmt16"

      elseif(isw.eq.1) then                ! Input material set data
        pstyp = ndm                        ! Sets plot dimension (1,2,3)

!       Read number pdes 
        nelu(:)   = 0
        du(:)     = 0
        dofu(:,:) = 0
        call elmt16_ISW01(v,npde,du,nelu,hist1,hist3)

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
        errck = pinput(td,3)
        d(1:3) = td(1:3)

!       Output material set data
        
        datades(1)="Em-Youngs modulus"

        datades(2)="[Nu]-Poisson ratio"

        datades(3)="[Rho]0-Density"

        write(iow,"(10x,f15.5,A3,A)")
     #     (d(i)," = ",datades(i),i=1,3)

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
          call elmt16_ISW03(v,d,xl,ua,k,c,m,r,hr(nh1),hr(nh2),gp,ngpo)

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
          call elmt16_ISW06(v,d,xl,ua,r,hr(nh1),hr(nh2),gp,ngpo)

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
        call elmt16_ISW05(v,d,xl,ua,m,r,hr(nh1),hr(nh2),gp,ngpo)

      elseif(isw.eq.9) then                ! Compute mass matrix

!       Map dofs
        call SB_ua_set(ul,ndf,nen,ua,du,dofu,nelu,npde)

!       Initialize element arrays
        c(:,:) = 0.0d0
        r(:)   = 0.0d0

!       Call the routine for ISW=9
        call elmt16_ISW09(v,d,xl,ua,c,r,hr(nh1),hr(nh2),gp,ngpo)

      endif

      end subroutine elmt16


!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt16_ISW01(v,npde,du,nelu,hist1,hist3)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER npde,du(10),nelu(10),hist1,hist3
      DOUBLE PRECISION v(1125)
      npde=(1)
      du(1)=(2)
      nelu(1)=(4)
      hist1=(0)
      hist3=(0)
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt16_ISW03(v,d,xl,ul,s,c,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i6
      DOUBLE PRECISION v(1125),d(3),xl(2,4),ul(6,8),s(8,8),c
     &(8,8),m(8,8),p(8),ht(0),hp(0),gp(4,4)
      v(73)=ul(1,8)
      v(72)=ul(1,7)
      v(71)=ul(1,6)
      v(70)=ul(1,5)
      v(69)=ul(1,4)
      v(68)=ul(1,3)
      v(67)=ul(1,2)
      v(66)=ul(1,1)
      v(17)=xl(2,4)
      v(16)=xl(1,4)
      v(15)=xl(2,3)
      v(984)=v(15)-v(17)
      v(14)=xl(1,3)
      v(986)=v(14)-v(16)
      v(13)=xl(2,2)
      v(980)=v(13)-v(15)
      v(12)=xl(1,2)
      v(982)=v(12)-v(14)
      v(11)=xl(2,1)
      v(983)=v(11)-v(13)
      v(979)=v(11)-v(17)
      v(10)=xl(1,1)
      v(985)=v(10)-v(12)
      v(981)=v(10)-v(16)
      v(9)=d(3)
      v(8)=d(2)
      v(965)=d(1)/(1d0+v(8))
      v(164)=(v(8)*v(965))/(1d0-2d0*v(8))
      v(162)=v(965)/2d0
      v(976)=2d0*v(162)
      v(987)=v(164)+v(976)
      DO i6=1,ngpo
       v(90)=gp(1,i6)
       v(99)=1d0-v(90)
       v(105)=(-0.25d0)*v(99)
       v(97)=1d0+v(90)
       v(106)=(-0.25d0)*v(97)
       v(111)=v(105)*v(979)+v(106)*v(980)
       v(109)=v(105)*v(981)+v(106)*v(982)
       v(91)=gp(2,i6)
       v(107)=(1d0+v(91))/4d0
       v(104)=((-1d0)+v(91))/4d0
       v(110)=v(104)*v(983)+v(107)*v(984)
       v(108)=v(104)*v(985)+v(107)*v(986)
       v(94)=-(v(104)*v(99))
       v(96)=-(v(104)*v(97))
       v(98)=v(107)*v(97)
       v(101)=v(107)*v(99)
       v(112)=-(v(109)*v(110))+v(108)*v(111)
       v(969)=gp(4,i6)*v(112)
       v(966)=v(9)*v(969)
       v(968)=v(966)*v(98)
       v(967)=v(94)*v(966)
       v(325)=(v(101)*v(101))*v(966)
       v(323)=v(966)*(v(98)*v(98))
       v(321)=v(966)*(v(96)*v(96))
       v(317)=(v(94)*v(94))*v(966)
       v(977)=v(9)*(ul(5,7)*v(101)+ul(5,1)*v(94)+ul(5,3)*v(96)+ul(5,5
     & )*v(98))
       v(978)=v(9)*(ul(5,8)*v(101)+ul(5,2)*v(94)+ul(5,4)*v(96)+ul(5,6
     & )*v(98))
       v(320)=v(101)*v(967)
       v(319)=v(967)*v(98)
       v(318)=v(96)*v(967)
       v(322)=v(96)*v(968)
       v(324)=v(101)*v(968)
       v(135)=-(v(111)/v(112))
       v(149)=-(v(107)*v(135))
       v(141)=-(v(104)*v(135))
       v(136)=v(109)/v(112)
       v(152)=-(v(107)*v(136))
       v(143)=-(v(104)*v(136))
       v(137)=-(v(110)/v(112))
       v(150)=v(105)*v(137)
       v(145)=v(106)*v(137)
       v(138)=v(108)/v(112)
       v(153)=v(105)*v(138)
       v(147)=v(106)*v(138)
       v(139)=v(141)+v(150)
       v(225)=v(139)*v(162)
       v(140)=v(143)+v(153)
       v(213)=v(140)*v(164)
       v(273)=(v(139)*v(213)+v(140)*v(225))*v(969)
       v(142)=-v(141)+v(145)
       v(970)=v(142)*v(162)
       v(214)=v(142)*v(164)
       v(221)=v(214)+2d0*v(970)
       v(144)=-v(143)+v(147)
       v(971)=v(144)*v(162)
       v(215)=v(144)*v(164)
       v(288)=v(969)*(v(142)*v(215)+v(144)*v(970))
       v(275)=v(969)*(v(139)*v(215)+v(140)*v(970))
       v(233)=v(215)+2d0*v(971)
       v(282)=v(969)*(v(140)*v(233)+v(139)*v(970))
       v(281)=v(969)*(v(140)*v(214)+v(139)*v(971))
       v(274)=v(969)*(v(139)*v(221)+v(140)*v(971))
       v(146)=-v(145)+v(149)
       v(972)=v(146)*v(162)
       v(216)=v(146)*v(164)
       v(222)=v(216)+2d0*v(972)
       v(148)=-v(147)+v(152)
       v(973)=v(148)*v(162)
       v(217)=v(148)*v(164)
       v(299)=v(969)*(v(146)*v(217)+v(148)*v(972))
       v(290)=v(969)*(v(142)*v(217)+v(144)*v(972))
       v(277)=v(969)*(v(139)*v(217)+v(140)*v(972))
       v(234)=v(217)+2d0*v(973)
       v(295)=v(969)*(v(144)*v(234)+v(142)*v(972))
       v(284)=v(969)*(v(140)*v(234)+v(139)*v(972))
       v(294)=v(969)*(v(144)*v(216)+v(142)*v(973))
       v(289)=v(969)*(v(142)*v(222)+v(144)*v(973))
       v(283)=v(969)*(v(140)*v(216)+v(139)*v(973))
       v(276)=v(969)*(v(139)*v(222)+v(140)*v(973))
       v(151)=-v(149)-v(150)
       v(974)=v(151)*v(162)
       v(218)=v(151)*v(164)
       v(223)=v(218)+2d0*v(974)
       v(154)=-v(152)-v(153)
       v(975)=v(154)*v(162)
       v(219)=v(154)*v(164)
       v(306)=v(969)*(v(151)*v(219)+v(154)*v(974))
       v(301)=v(969)*(v(146)*v(219)+v(148)*v(974))
       v(292)=v(969)*(v(142)*v(219)+v(144)*v(974))
       v(279)=v(969)*(v(139)*v(219)+v(140)*v(974))
       v(235)=v(219)+2d0*v(975)
       v(304)=v(969)*(v(148)*v(235)+v(146)*v(974))
       v(297)=v(969)*(v(144)*v(235)+v(142)*v(974))
       v(286)=v(969)*(v(140)*v(235)+v(139)*v(974))
       v(303)=v(969)*(v(148)*v(218)+v(146)*v(975))
       v(300)=v(969)*(v(146)*v(223)+v(148)*v(975))
       v(296)=v(969)*(v(144)*v(218)+v(142)*v(975))
       v(291)=v(969)*(v(142)*v(223)+v(144)*v(975))
       v(285)=v(969)*(v(140)*v(218)+v(139)*v(975))
       v(278)=v(969)*(v(139)*v(223)+v(140)*v(975))
       v(159)=v(139)*v(66)+v(142)*v(68)+v(146)*v(70)+v(151)*v(72)
       v(161)=v(140)*v(67)+v(144)*v(69)+v(148)*v(71)+v(154)*v(73)
       v(166)=(v(159)+v(161))*v(164)
       v(167)=v(166)+v(159)*v(976)
       v(168)=v(162)*(v(140)*v(66)+v(139)*v(67)+v(144)*v(68)+v(142)*v
     & (69)+v(148)*v(70)+v(146)*v(71)+v(154)*v(72)+v(151)*v(73))
       v(169)=v(166)+v(161)*v(976)
       p(1)=p(1)+v(969)*(v(139)*v(167)+v(140)*v(168)+v(94)*v(977))
       p(2)=p(2)+v(969)*(v(139)*v(168)+v(140)*v(169)+v(94)*v(978))
       p(3)=p(3)+v(969)*(v(142)*v(167)+v(144)*v(168)+v(96)*v(977))
       p(4)=p(4)+v(969)*(v(142)*v(168)+v(144)*v(169)+v(96)*v(978))
       p(5)=p(5)+v(969)*(v(146)*v(167)+v(148)*v(168)+v(977)*v(98))
       p(6)=p(6)+v(969)*(v(146)*v(168)+v(148)*v(169)+v(978)*v(98))
       p(7)=p(7)+v(969)*(v(151)*v(167)+v(154)*v(168)+v(101)*v(977))
       p(8)=p(8)+v(969)*(v(151)*v(168)+v(154)*v(169)+v(101)*v(978))
       s(1,1)=s(1,1)+v(969)*((v(140)*v(140))*v(162)+(v(139)*v(139))*v
     & (987))
       s(1,2)=s(1,2)+v(273)
       s(1,3)=s(1,3)+v(274)
       s(1,4)=s(1,4)+v(275)
       s(1,5)=s(1,5)+v(276)
       s(1,6)=s(1,6)+v(277)
       s(1,7)=s(1,7)+v(278)
       s(1,8)=s(1,8)+v(279)
       s(2,1)=s(2,1)+v(273)
       s(2,2)=s(2,2)+v(969)*(v(139)*v(225)+v(140)*(v(213)+v(140)*v
     & (976)))
       s(2,3)=s(2,3)+v(281)
       s(2,4)=s(2,4)+v(282)
       s(2,5)=s(2,5)+v(283)
       s(2,6)=s(2,6)+v(284)
       s(2,7)=s(2,7)+v(285)
       s(2,8)=s(2,8)+v(286)
       s(3,1)=s(3,1)+v(274)
       s(3,2)=s(3,2)+v(281)
       s(3,3)=s(3,3)+v(969)*(v(142)*v(221)+v(144)*v(971))
       s(3,4)=s(3,4)+v(288)
       s(3,5)=s(3,5)+v(289)
       s(3,6)=s(3,6)+v(290)
       s(3,7)=s(3,7)+v(291)
       s(3,8)=s(3,8)+v(292)
       s(4,1)=s(4,1)+v(275)
       s(4,2)=s(4,2)+v(282)
       s(4,3)=s(4,3)+v(288)
       s(4,4)=s(4,4)+v(969)*(v(144)*v(233)+v(142)*v(970))
       s(4,5)=s(4,5)+v(294)
       s(4,6)=s(4,6)+v(295)
       s(4,7)=s(4,7)+v(296)
       s(4,8)=s(4,8)+v(297)
       s(5,1)=s(5,1)+v(276)
       s(5,2)=s(5,2)+v(283)
       s(5,3)=s(5,3)+v(289)
       s(5,4)=s(5,4)+v(294)
       s(5,5)=s(5,5)+v(969)*(v(146)*v(222)+v(148)*v(973))
       s(5,6)=s(5,6)+v(299)
       s(5,7)=s(5,7)+v(300)
       s(5,8)=s(5,8)+v(301)
       s(6,1)=s(6,1)+v(277)
       s(6,2)=s(6,2)+v(284)
       s(6,3)=s(6,3)+v(290)
       s(6,4)=s(6,4)+v(295)
       s(6,5)=s(6,5)+v(299)
       s(6,6)=s(6,6)+v(969)*(v(148)*v(234)+v(146)*v(972))
       s(6,7)=s(6,7)+v(303)
       s(6,8)=s(6,8)+v(304)
       s(7,1)=s(7,1)+v(278)
       s(7,2)=s(7,2)+v(285)
       s(7,3)=s(7,3)+v(291)
       s(7,4)=s(7,4)+v(296)
       s(7,5)=s(7,5)+v(300)
       s(7,6)=s(7,6)+v(303)
       s(7,7)=s(7,7)+v(969)*(v(151)*v(223)+v(154)*v(975))
       s(7,8)=s(7,8)+v(306)
       s(8,1)=s(8,1)+v(279)
       s(8,2)=s(8,2)+v(286)
       s(8,3)=s(8,3)+v(292)
       s(8,4)=s(8,4)+v(297)
       s(8,5)=s(8,5)+v(301)
       s(8,6)=s(8,6)+v(304)
       s(8,7)=s(8,7)+v(306)
       s(8,8)=s(8,8)+v(969)*(v(154)*v(235)+v(151)*v(974))
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
       m(1,1)=m(1,1)+v(317)
       m(1,2)=0d0+m(1,2)
       m(1,3)=m(1,3)+v(318)
       m(1,4)=0d0+m(1,4)
       m(1,5)=m(1,5)+v(319)
       m(1,6)=0d0+m(1,6)
       m(1,7)=m(1,7)+v(320)
       m(1,8)=0d0+m(1,8)
       m(2,1)=0d0+m(2,1)
       m(2,2)=m(2,2)+v(317)
       m(2,3)=0d0+m(2,3)
       m(2,4)=m(2,4)+v(318)
       m(2,5)=0d0+m(2,5)
       m(2,6)=m(2,6)+v(319)
       m(2,7)=0d0+m(2,7)
       m(2,8)=m(2,8)+v(320)
       m(3,1)=m(3,1)+v(318)
       m(3,2)=0d0+m(3,2)
       m(3,3)=m(3,3)+v(321)
       m(3,4)=0d0+m(3,4)
       m(3,5)=m(3,5)+v(322)
       m(3,6)=0d0+m(3,6)
       m(3,7)=m(3,7)+v(319)
       m(3,8)=0d0+m(3,8)
       m(4,1)=0d0+m(4,1)
       m(4,2)=m(4,2)+v(318)
       m(4,3)=0d0+m(4,3)
       m(4,4)=m(4,4)+v(321)
       m(4,5)=0d0+m(4,5)
       m(4,6)=m(4,6)+v(322)
       m(4,7)=0d0+m(4,7)
       m(4,8)=m(4,8)+v(319)
       m(5,1)=m(5,1)+v(319)
       m(5,2)=0d0+m(5,2)
       m(5,3)=m(5,3)+v(322)
       m(5,4)=0d0+m(5,4)
       m(5,5)=m(5,5)+v(323)
       m(5,6)=0d0+m(5,6)
       m(5,7)=m(5,7)+v(324)
       m(5,8)=0d0+m(5,8)
       m(6,1)=0d0+m(6,1)
       m(6,2)=m(6,2)+v(319)
       m(6,3)=0d0+m(6,3)
       m(6,4)=m(6,4)+v(322)
       m(6,5)=0d0+m(6,5)
       m(6,6)=m(6,6)+v(323)
       m(6,7)=0d0+m(6,7)
       m(6,8)=m(6,8)+v(324)
       m(7,1)=m(7,1)+v(320)
       m(7,2)=0d0+m(7,2)
       m(7,3)=m(7,3)+v(319)
       m(7,4)=0d0+m(7,4)
       m(7,5)=m(7,5)+v(324)
       m(7,6)=0d0+m(7,6)
       m(7,7)=m(7,7)+v(325)
       m(7,8)=0d0+m(7,8)
       m(8,1)=0d0+m(8,1)
       m(8,2)=m(8,2)+v(320)
       m(8,3)=0d0+m(8,3)
       m(8,4)=m(8,4)+v(319)
       m(8,5)=0d0+m(8,5)
       m(8,6)=m(8,6)+v(324)
       m(8,7)=0d0+m(8,7)
       m(8,8)=m(8,8)+v(325)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt16_ISW05(v,d,xl,ul,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i329
      DOUBLE PRECISION v(1125),d(3),xl(2,4),ul(6,8),m(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      v(340)=xl(2,4)
      v(339)=xl(1,4)
      v(338)=xl(2,3)
      v(998)=v(338)-v(340)
      v(337)=xl(1,3)
      v(994)=v(337)-v(339)
      v(336)=xl(2,2)
      v(992)=v(336)-v(338)
      v(335)=xl(1,2)
      v(996)=v(335)-v(337)
      v(334)=xl(2,1)
      v(997)=v(334)-v(336)
      v(991)=v(334)-v(340)
      v(333)=xl(1,1)
      v(995)=v(333)-v(339)
      v(993)=v(333)-v(335)
      DO i329=1,ngpo
       v(413)=gp(1,i329)
       v(422)=1d0-v(413)
       v(428)=(-0.25d0)*v(422)
       v(420)=1d0+v(413)
       v(429)=(-0.25d0)*v(420)
       v(414)=gp(2,i329)
       v(430)=(1d0+v(414))/4d0
       v(427)=((-1d0)+v(414))/4d0
       v(417)=-(v(422)*v(427))
       v(419)=-(v(420)*v(427))
       v(421)=v(420)*v(430)
       v(424)=v(422)*v(430)
       v(988)=d(3)*gp(4,i329)*((v(428)*v(991)+v(429)*v(992))*(v(427
     & )*v(993)+v(430)*v(994))-(v(428)*v(995)+v(429)*v(996))*(v(427
     & )*v(997)+v(430)*v(998)))
       v(990)=v(421)*v(988)
       v(989)=v(417)*v(988)
       v(551)=(v(424)*v(424))*v(988)
       v(549)=(v(421)*v(421))*v(988)
       v(547)=(v(419)*v(419))*v(988)
       v(543)=(v(417)*v(417))*v(988)
       v(546)=v(424)*v(989)
       v(545)=v(421)*v(989)
       v(544)=v(419)*v(989)
       v(548)=v(419)*v(990)
       v(550)=v(424)*v(990)
       m(1,1)=m(1,1)+v(543)
       m(1,2)=0d0+m(1,2)
       m(1,3)=m(1,3)+v(544)
       m(1,4)=0d0+m(1,4)
       m(1,5)=m(1,5)+v(545)
       m(1,6)=0d0+m(1,6)
       m(1,7)=m(1,7)+v(546)
       m(1,8)=0d0+m(1,8)
       m(2,1)=0d0+m(2,1)
       m(2,2)=m(2,2)+v(543)
       m(2,3)=0d0+m(2,3)
       m(2,4)=m(2,4)+v(544)
       m(2,5)=0d0+m(2,5)
       m(2,6)=m(2,6)+v(545)
       m(2,7)=0d0+m(2,7)
       m(2,8)=m(2,8)+v(546)
       m(3,1)=m(3,1)+v(544)
       m(3,2)=0d0+m(3,2)
       m(3,3)=m(3,3)+v(547)
       m(3,4)=0d0+m(3,4)
       m(3,5)=m(3,5)+v(548)
       m(3,6)=0d0+m(3,6)
       m(3,7)=m(3,7)+v(545)
       m(3,8)=0d0+m(3,8)
       m(4,1)=0d0+m(4,1)
       m(4,2)=m(4,2)+v(544)
       m(4,3)=0d0+m(4,3)
       m(4,4)=m(4,4)+v(547)
       m(4,5)=0d0+m(4,5)
       m(4,6)=m(4,6)+v(548)
       m(4,7)=0d0+m(4,7)
       m(4,8)=m(4,8)+v(545)
       m(5,1)=m(5,1)+v(545)
       m(5,2)=0d0+m(5,2)
       m(5,3)=m(5,3)+v(548)
       m(5,4)=0d0+m(5,4)
       m(5,5)=m(5,5)+v(549)
       m(5,6)=0d0+m(5,6)
       m(5,7)=m(5,7)+v(550)
       m(5,8)=0d0+m(5,8)
       m(6,1)=0d0+m(6,1)
       m(6,2)=m(6,2)+v(545)
       m(6,3)=0d0+m(6,3)
       m(6,4)=m(6,4)+v(548)
       m(6,5)=0d0+m(6,5)
       m(6,6)=m(6,6)+v(549)
       m(6,7)=0d0+m(6,7)
       m(6,8)=m(6,8)+v(550)
       m(7,1)=m(7,1)+v(546)
       m(7,2)=0d0+m(7,2)
       m(7,3)=m(7,3)+v(545)
       m(7,4)=0d0+m(7,4)
       m(7,5)=m(7,5)+v(550)
       m(7,6)=0d0+m(7,6)
       m(7,7)=m(7,7)+v(551)
       m(7,8)=0d0+m(7,8)
       m(8,1)=0d0+m(8,1)
       m(8,2)=m(8,2)+v(546)
       m(8,3)=0d0+m(8,3)
       m(8,4)=m(8,4)+v(545)
       m(8,5)=0d0+m(8,5)
       m(8,6)=m(8,6)+v(550)
       m(8,7)=0d0+m(8,7)
       m(8,8)=m(8,8)+v(551)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt16_ISW06(v,d,xl,ul,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i553
      DOUBLE PRECISION v(1125),d(3),xl(2,4),ul(6,8),p(8),ht
     &(0),hp(0),gp(4,4)
      v(620)=ul(1,8)
      v(619)=ul(1,7)
      v(618)=ul(1,6)
      v(617)=ul(1,5)
      v(616)=ul(1,4)
      v(615)=ul(1,3)
      v(614)=ul(1,2)
      v(613)=ul(1,1)
      v(564)=xl(2,4)
      v(563)=xl(1,4)
      v(562)=xl(2,3)
      v(1009)=v(562)-v(564)
      v(561)=xl(1,3)
      v(1011)=v(561)-v(563)
      v(560)=xl(2,2)
      v(1005)=v(560)-v(562)
      v(559)=xl(1,2)
      v(1007)=v(559)-v(561)
      v(558)=xl(2,1)
      v(1008)=v(558)-v(560)
      v(1004)=v(558)-v(564)
      v(557)=xl(1,1)
      v(1010)=v(557)-v(559)
      v(1006)=v(557)-v(563)
      v(556)=d(3)
      v(555)=d(2)
      v(999)=d(1)/(1d0+v(555))
      v(711)=(v(555)*v(999))/(1d0-2d0*v(555))
      v(709)=v(999)/2d0
      v(1000)=2d0*v(709)
      DO i553=1,ngpo
       v(637)=gp(1,i553)
       v(646)=1d0-v(637)
       v(652)=(-0.25d0)*v(646)
       v(644)=1d0+v(637)
       v(653)=(-0.25d0)*v(644)
       v(658)=v(1004)*v(652)+v(1005)*v(653)
       v(656)=v(1006)*v(652)+v(1007)*v(653)
       v(638)=gp(2,i553)
       v(654)=(1d0+v(638))/4d0
       v(651)=((-1d0)+v(638))/4d0
       v(657)=v(1008)*v(651)+v(1009)*v(654)
       v(655)=v(1010)*v(651)+v(1011)*v(654)
       v(641)=-(v(646)*v(651))
       v(643)=-(v(644)*v(651))
       v(645)=v(644)*v(654)
       v(648)=v(646)*v(654)
       v(659)=-(v(656)*v(657))+v(655)*v(658)
       v(1001)=gp(4,i553)*v(659)
       v(1002)=v(556)*(ul(5,1)*v(641)+ul(5,3)*v(643)+ul(5,5)*v(645)
     & +ul(5,7)*v(648))
       v(1003)=v(556)*(ul(5,2)*v(641)+ul(5,4)*v(643)+ul(5,6)*v(645)
     & +ul(5,8)*v(648))
       v(682)=-(v(658)/v(659))
       v(696)=-(v(654)*v(682))
       v(688)=-(v(651)*v(682))
       v(683)=v(656)/v(659)
       v(699)=-(v(654)*v(683))
       v(690)=-(v(651)*v(683))
       v(684)=-(v(657)/v(659))
       v(697)=v(652)*v(684)
       v(692)=v(653)*v(684)
       v(685)=v(655)/v(659)
       v(700)=v(652)*v(685)
       v(694)=v(653)*v(685)
       v(686)=v(688)+v(697)
       v(687)=v(690)+v(700)
       v(689)=-v(688)+v(692)
       v(691)=-v(690)+v(694)
       v(693)=-v(692)+v(696)
       v(695)=-v(694)+v(699)
       v(698)=-v(696)-v(697)
       v(701)=-v(699)-v(700)
       v(706)=v(613)*v(686)+v(615)*v(689)+v(617)*v(693)+v(619)*v(698)
       v(708)=v(614)*v(687)+v(616)*v(691)+v(618)*v(695)+v(620)*v(701)
       v(713)=(v(706)+v(708))*v(711)
       v(714)=v(1000)*v(706)+v(713)
       v(715)=(v(614)*v(686)+v(613)*v(687)+v(616)*v(689)+v(615)*v(691
     & )+v(618)*v(693)+v(617)*v(695)+v(620)*v(698)+v(619)*v(701))*v
     & (709)
       v(716)=v(1000)*v(708)+v(713)
       p(1)=p(1)+v(1001)*(v(1002)*v(641)+v(686)*v(714)+v(687)*v(715))
       p(2)=p(2)+v(1001)*(v(1003)*v(641)+v(686)*v(715)+v(687)*v(716))
       p(3)=p(3)+v(1001)*(v(1002)*v(643)+v(689)*v(714)+v(691)*v(715))
       p(4)=p(4)+v(1001)*(v(1003)*v(643)+v(689)*v(715)+v(691)*v(716))
       p(5)=p(5)+v(1001)*(v(1002)*v(645)+v(693)*v(714)+v(695)*v(715))
       p(6)=p(6)+v(1001)*(v(1003)*v(645)+v(693)*v(715)+v(695)*v(716))
       p(7)=p(7)+v(1001)*(v(1002)*v(648)+v(698)*v(714)+v(701)*v(715))
       p(8)=p(8)+v(1001)*(v(1003)*v(648)+v(698)*v(715)+v(701)*v(716))
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt16_ISW09(v,d,xl,ul,c,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i759
      DOUBLE PRECISION v(1125),d(3),xl(2,4),ul(6,8),c(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i759=1,ngpo
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
