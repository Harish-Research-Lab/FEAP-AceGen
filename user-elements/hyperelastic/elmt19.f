!**************************************************************
!* AceGen    7.205 MacOSX (15 Jan 21)                         *
!*           Co. J. Korelc  2020           8 Aug 21 00:11:17  *
!**************************************************************
! User     : Full professional version
! Notebook : Yeoh
! Evaluation time                 : 10 s    Mode  : Optimal
! Number of formulae              : 281     Method: Automatic
! Subroutine                      : elmt19_ISW01 size: 57
! Subroutine                      : elmt19_ISW03 size: 4634
! Subroutine                      : elmt19_ISW05 size: 450
! Subroutine                      : elmt19_ISW06 size: 1016
! Subroutine                      : elmt19_ISW09 size: 450
! Total size of Mathematica  code : 6607 subexpressions
! Total size of Fortran code      : 21716 bytes

      subroutine elmt19 (d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

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

      character (len=50) :: datades(4)
      logical       ::  pinput
      integer       ::  dofacegen,ngpo
      real (kind=8) ::  td(10)
      real (kind=8) ::  v(1244)   ! AceGen Storage
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

      symmetric=.true.

      if(isw.ne.1) then

!       Read number of PDEs
        npde = nint(d(100))

      endif
 
      if(isw.lt.0) then
        utx(1) = "elmt19"

      elseif(isw.eq.1) then                ! Input material set data
        pstyp = ndm                        ! Sets plot dimension (1,2,3)

!       Read number pdes 
        nelu(:)   = 0
        du(:)     = 0
        dofu(:,:) = 0
        call elmt19_ISW01(v,npde,du,nelu,hist1,hist3)

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
        errck = pinput(td,4)
        d(1:4) = td(1:4)

!       Output material set data
        
        datades(1)="Em-Youngs modulus"

        datades(2)="[Nu]-Poisson ratio"

        datades(3)="k1p-Para01"

        datades(4)="k2p-Para02"

        write(iow,"(10x,f15.5,A3,A)")
     #     (d(i)," = ",datades(i),i=1,4)

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
          call elmt19_ISW03(v,d,xl,ua,k,c,m,r,hr(nh1),hr(nh2),gp,ngpo)

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
          call elmt19_ISW06(v,d,xl,ua,r,hr(nh1),hr(nh2),gp,ngpo)

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
        call elmt19_ISW05(v,d,xl,ua,m,r,hr(nh1),hr(nh2),gp,ngpo)

      elseif(isw.eq.9) then                ! Compute damping matrix

!       Map dofs
        call SB_ua_set(ul,ndf,nen,ua,du,dofu,nelu,npde)

!       Initialize element arrays
        c(:,:) = 0.0d0
        r(:)   = 0.0d0

!       Call the routine for ISW=9
        call elmt19_ISW09(v,d,xl,ua,c,r,hr(nh1),hr(nh2),gp,ngpo)

      endif

      end subroutine elmt19


!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt19_ISW01(v,npde,du,nelu,hist1,hist3)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER npde,du(10),nelu(10),hist1,hist3
      DOUBLE PRECISION v(1244)
      npde=(1)
      du(1)=(2)
      nelu(1)=(4)
      hist1=(0)
      hist3=(0)
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt19_ISW03(v,d,xl,ul,s,c,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i6
      DOUBLE PRECISION v(1244),d(2),xl(2,4),ul(6,8),s(8,8),c
     &(8,8),m(8,8),p(8),ht(0),hp(0),gp(4,4)
      v(74)=ul(1,8)
      v(73)=ul(1,7)
      v(72)=ul(1,6)
      v(71)=ul(1,5)
      v(70)=ul(1,4)
      v(69)=ul(1,3)
      v(68)=ul(1,2)
      v(67)=ul(1,1)
      v(18)=xl(2,4)
      v(17)=xl(1,4)
      v(16)=xl(2,3)
      v(1119)=v(16)-v(18)
      v(15)=xl(1,3)
      v(1121)=v(15)-v(17)
      v(14)=xl(2,2)
      v(1115)=v(14)-v(16)
      v(13)=xl(1,2)
      v(1117)=v(13)-v(15)
      v(12)=xl(2,1)
      v(1118)=v(12)-v(14)
      v(1114)=v(12)-v(18)
      v(11)=xl(1,1)
      v(1120)=v(11)-v(13)
      v(1116)=v(11)-v(17)
      v(1111)=3d0*d(4)
      v(9)=d(3)
      v(8)=d(2)
      v(7)=d(1)
      v(188)=v(7)/(2d0+2d0*v(8))
      v(187)=v(7)/(3d0-6d0*v(8))
      DO i6=1,ngpo
       v(91)=gp(1,i6)
       v(106)=((-1d0)+v(91))/4d0
       v(107)=((-1d0)-v(91))/4d0
       v(112)=v(106)*v(1114)+v(107)*v(1115)
       v(110)=v(106)*v(1116)+v(107)*v(1117)
       v(92)=gp(2,i6)
       v(108)=(1d0+v(92))/4d0
       v(105)=((-1d0)+v(92))/4d0
       v(111)=v(105)*v(1118)+v(108)*v(1119)
       v(109)=v(105)*v(1120)+v(108)*v(1121)
       v(113)=-(v(110)*v(111))+v(109)*v(112)
       v(1113)=gp(4,i6)*v(113)
       v(136)=-(v(112)/v(113))
       v(150)=-(v(108)*v(136))
       v(142)=-(v(105)*v(136))
       v(137)=v(110)/v(113)
       v(153)=-(v(108)*v(137))
       v(144)=-(v(105)*v(137))
       v(138)=-(v(111)/v(113))
       v(151)=v(106)*v(138)
       v(146)=v(107)*v(138)
       v(139)=v(109)/v(113)
       v(154)=v(106)*v(139)
       v(148)=v(107)*v(139)
       v(140)=v(142)+v(151)
       v(141)=v(144)+v(154)
       v(143)=-v(142)+v(146)
       v(145)=-v(144)+v(148)
       v(147)=-v(146)+v(150)
       v(149)=-v(148)+v(153)
       v(152)=-v(150)-v(151)
       v(155)=-v(153)-v(154)
       v(164)=1d0+v(140)*v(67)+v(143)*v(69)+v(147)*v(71)+v(152)*v(73)
       v(165)=v(141)*v(67)+v(145)*v(69)+v(149)*v(71)+v(155)*v(73)
       v(260)=v(155)*v(164)-v(152)*v(165)
       v(258)=v(149)*v(164)-v(147)*v(165)
       v(256)=v(145)*v(164)-v(143)*v(165)
       v(254)=v(141)*v(164)-v(140)*v(165)
       v(251)=2d0*(v(152)*v(164)+v(155)*v(165))
       v(249)=2d0*(v(147)*v(164)+v(149)*v(165))
       v(247)=2d0*(v(143)*v(164)+v(145)*v(165))
       v(245)=2d0*(v(140)*v(164)+v(141)*v(165))
       v(167)=v(140)*v(68)+v(143)*v(70)+v(147)*v(72)+v(152)*v(74)
       v(168)=1d0+v(141)*v(68)+v(145)*v(70)+v(149)*v(72)+v(155)*v(74)
       v(259)=-(v(155)*v(167))+v(152)*v(168)
       v(257)=-(v(149)*v(167))+v(147)*v(168)
       v(255)=-(v(145)*v(167))+v(143)*v(168)
       v(253)=-(v(141)*v(167))+v(140)*v(168)
       v(252)=2d0*(v(152)*v(167)+v(155)*v(168))
       v(250)=2d0*(v(147)*v(167)+v(149)*v(168))
       v(248)=2d0*(v(143)*v(167)+v(145)*v(168))
       v(246)=2d0*(v(140)*v(167)+v(141)*v(168))
       v(199)=1d0+(v(164)*v(164))+(v(165)*v(165))+(v(167)*v(167))+(v
     & (168)*v(168))
       v(185)=-(v(165)*v(167))+v(164)*v(168)
       v(1109)=(-2d0/3d0)/v(185)**0.16666666666666669d1
       v(305)=v(1109)*v(199)
       v(269)=v(1109)*v(260)
       v(268)=v(1109)*v(259)
       v(267)=v(1109)*v(258)
       v(266)=v(1109)*v(257)
       v(265)=v(1109)*v(256)
       v(264)=v(1109)*v(255)
       v(263)=v(1109)*v(254)
       v(262)=v(1109)*v(253)
       v(201)=1d0/v(185)**0.6666666666666666d0
       v(193)=(-3d0)+v(199)*v(201)
       v(1110)=v(188)*(v(1111)*v(193)+v(9))
       v(293)=v(1110)*(v(201)*v(252)+v(199)*v(269))
       v(292)=v(1110)*(v(201)*v(251)+v(199)*v(268))
       v(291)=v(1110)*(v(201)*v(250)+v(199)*v(267))
       v(290)=v(1110)*(v(201)*v(249)+v(199)*v(266))
       v(289)=v(1110)*(v(201)*v(248)+v(199)*v(265))
       v(288)=v(1110)*(v(201)*v(247)+v(199)*v(264))
       v(287)=v(1110)*(v(201)*v(246)+v(199)*v(263))
       v(286)=v(1110)*(v(201)*v(245)+v(199)*v(262))
       v(197)=(v(188)*(1d0+v(1111)*(v(193)*v(193))+2d0*v(193)*v(9)))
     & /2d0
       v(1112)=v(187)+((10d0/9d0)*v(197)*v(199))/v(185
     & )**0.26666666666666666d1
       v(302)=v(1109)*v(197)
       v(312)=v(1112)*v(260)+v(252)*v(302)+v(293)*v(305)
       v(311)=v(1112)*v(259)+v(251)*v(302)+v(292)*v(305)
       v(310)=v(1112)*v(258)+v(250)*v(302)+v(291)*v(305)
       v(309)=v(1112)*v(257)+v(249)*v(302)+v(290)*v(305)
       v(308)=v(1112)*v(256)+v(248)*v(302)+v(289)*v(305)
       v(307)=v(1112)*v(255)+v(247)*v(302)+v(288)*v(305)
       v(306)=v(1112)*v(254)+v(246)*v(302)+v(287)*v(305)
       v(303)=v(1112)*v(253)+v(245)*v(302)+v(286)*v(305)
       v(301)=2d0*(v(197)*v(269)+v(201)*v(293))
       v(300)=2d0*(v(197)*v(268)+v(201)*v(292))
       v(299)=2d0*(v(197)*v(267)+v(201)*v(291))
       v(298)=2d0*(v(197)*v(266)+v(201)*v(290))
       v(297)=2d0*(v(197)*v(265)+v(201)*v(289))
       v(296)=2d0*(v(197)*v(264)+v(201)*v(288))
       v(295)=2d0*(v(197)*v(263)+v(201)*v(287))
       v(294)=2d0*(v(197)*v(262)+v(201)*v(286))
       v(203)=2d0*v(197)*v(201)
       v(363)=v(155)*v(203)
       v(364)=v(168)*v(301)+v(164)*v(312)+v(363)
       v(359)=v(149)*v(203)
       v(360)=v(168)*v(299)+v(164)*v(310)+v(359)
       v(355)=v(145)*v(203)
       v(356)=v(168)*v(297)+v(164)*v(308)+v(355)
       v(351)=v(141)*v(203)
       v(346)=v(152)*v(203)
       v(347)=v(167)*v(301)-v(165)*v(312)+v(346)
       v(433)=v(1113)*(v(147)*v(347)+v(149)*v(364))
       v(426)=v(1113)*(v(143)*v(347)+v(145)*v(364))
       v(415)=v(1113)*(v(140)*v(347)+v(141)*v(364))
       v(342)=v(147)*v(203)
       v(343)=v(167)*v(299)-v(165)*v(310)+v(342)
       v(424)=v(1113)*(v(143)*v(343)+v(145)*v(360))
       v(413)=v(1113)*(v(140)*v(343)+v(141)*v(360))
       v(338)=v(143)*v(203)
       v(339)=v(167)*v(297)-v(165)*v(308)+v(338)
       v(411)=v(1113)*(v(140)*v(339)+v(141)*v(356))
       v(334)=v(140)*v(203)
       v(329)=v(165)*v(300)-v(167)*v(311)+v(363)
       v(327)=v(165)*v(298)-v(167)*v(309)+v(359)
       v(325)=v(165)*v(296)-v(167)*v(307)+v(355)
       v(320)=v(164)*v(300)+v(168)*v(311)+v(346)
       v(429)=v(1113)*(v(147)*v(320)+v(149)*v(329))
       v(420)=v(1113)*(v(143)*v(320)+v(145)*v(329))
       v(407)=v(1113)*(v(140)*v(320)+v(141)*v(329))
       v(318)=v(164)*v(298)+v(168)*v(309)+v(342)
       v(418)=v(1113)*(v(143)*v(318)+v(145)*v(327))
       v(405)=v(1113)*(v(140)*v(318)+v(141)*v(327))
       v(316)=v(164)*v(296)+v(168)*v(307)+v(338)
       v(403)=v(1113)*(v(140)*v(316)+v(141)*v(325))
       v(200)=((-1d0)+v(185))*v(187)+v(199)*v(302)
       v(361)=-(v(152)*v(200))
       v(362)=v(168)*v(300)+v(164)*v(311)-v(361)
       v(357)=-(v(147)*v(200))
       v(358)=v(168)*v(298)+v(164)*v(309)-v(357)
       v(353)=-(v(143)*v(200))
       v(344)=v(155)*v(200)
       v(345)=v(167)*v(300)-v(165)*v(311)-v(344)
       v(432)=v(1113)*(v(147)*v(345)+v(149)*v(362))
       v(425)=v(1113)*(v(143)*v(345)+v(145)*v(362))
       v(414)=v(1113)*(v(140)*v(345)+v(141)*v(362))
       v(340)=v(149)*v(200)
       v(341)=v(167)*v(298)-v(165)*v(309)-v(340)
       v(423)=v(1113)*(v(143)*v(341)+v(145)*v(358))
       v(412)=v(1113)*(v(140)*v(341)+v(141)*v(358))
       v(336)=v(145)*v(200)
       v(410)=v(1113)*(v(140)*(v(167)*v(296)-v(165)*v(307)-v(336))+v
     & (141)*(v(168)*v(296)+v(164)*v(307)-v(353)))
       v(330)=v(165)*v(301)-v(167)*v(312)+v(361)
       v(328)=v(165)*v(299)-v(167)*v(310)+v(357)
       v(326)=v(165)*v(297)-v(167)*v(308)+v(353)
       v(321)=v(164)*v(301)+v(168)*v(312)+v(344)
       v(435)=v(1113)*(v(152)*v(321)+v(155)*v(330))
       v(430)=v(1113)*(v(147)*v(321)+v(149)*v(330))
       v(421)=v(1113)*(v(143)*v(321)+v(145)*v(330))
       v(408)=v(1113)*(v(140)*v(321)+v(141)*v(330))
       v(319)=v(164)*v(299)+v(168)*v(310)+v(340)
       v(428)=v(1113)*(v(147)*v(319)+v(149)*v(328))
       v(419)=v(1113)*(v(143)*v(319)+v(145)*v(328))
       v(406)=v(1113)*(v(140)*v(319)+v(141)*v(328))
       v(317)=v(164)*v(297)+v(168)*v(308)+v(336)
       v(417)=v(1113)*(v(143)*v(317)+v(145)*v(326))
       v(404)=v(1113)*(v(140)*v(317)+v(141)*v(326))
       v(402)=v(1113)*(v(141)*(v(165)*v(295)-v(167)*v(306))+v(140)*(v
     & (164)*v(295)+v(168)*v(306)))
       v(202)=v(168)*v(200)+v(164)*v(203)
       v(204)=-(v(167)*v(200))+v(165)*v(203)
       v(206)=-(v(165)*v(200))+v(167)*v(203)
       v(207)=v(164)*v(200)+v(168)*v(203)
       p(1)=p(1)+v(1113)*(v(140)*v(202)+v(141)*v(204))
       p(2)=p(2)+v(1113)*(v(140)*v(206)+v(141)*v(207))
       p(3)=p(3)+v(1113)*(v(143)*v(202)+v(145)*v(204))
       p(4)=p(4)+v(1113)*(v(143)*v(206)+v(145)*v(207))
       p(5)=p(5)+v(1113)*(v(147)*v(202)+v(149)*v(204))
       p(6)=p(6)+v(1113)*(v(147)*v(206)+v(149)*v(207))
       p(7)=p(7)+v(1113)*(v(152)*v(202)+v(155)*v(204))
       p(8)=p(8)+v(1113)*(v(152)*v(206)+v(155)*v(207))
       s(1,1)=s(1,1)+v(1113)*(v(140)*(v(164)*v(294)+v(168)*v(303)+v
     & (334))+v(141)*(v(165)*v(294)-v(167)*v(303)+v(351)))
       s(1,2)=s(1,2)+v(402)
       s(1,3)=s(1,3)+v(403)
       s(1,4)=s(1,4)+v(404)
       s(1,5)=s(1,5)+v(405)
       s(1,6)=s(1,6)+v(406)
       s(1,7)=s(1,7)+v(407)
       s(1,8)=s(1,8)+v(408)
       s(2,1)=s(2,1)+v(402)
       s(2,2)=s(2,2)+v(1113)*(v(140)*(v(167)*v(295)-v(165)*v(306)+v
     & (334))+v(141)*(v(168)*v(295)+v(164)*v(306)+v(351)))
       s(2,3)=s(2,3)+v(410)
       s(2,4)=s(2,4)+v(411)
       s(2,5)=s(2,5)+v(412)
       s(2,6)=s(2,6)+v(413)
       s(2,7)=s(2,7)+v(414)
       s(2,8)=s(2,8)+v(415)
       s(3,1)=s(3,1)+v(403)
       s(3,2)=s(3,2)+v(410)
       s(3,3)=s(3,3)+v(1113)*(v(143)*v(316)+v(145)*v(325))
       s(3,4)=s(3,4)+v(417)
       s(3,5)=s(3,5)+v(418)
       s(3,6)=s(3,6)+v(419)
       s(3,7)=s(3,7)+v(420)
       s(3,8)=s(3,8)+v(421)
       s(4,1)=s(4,1)+v(404)
       s(4,2)=s(4,2)+v(411)
       s(4,3)=s(4,3)+v(417)
       s(4,4)=s(4,4)+v(1113)*(v(143)*v(339)+v(145)*v(356))
       s(4,5)=s(4,5)+v(423)
       s(4,6)=s(4,6)+v(424)
       s(4,7)=s(4,7)+v(425)
       s(4,8)=s(4,8)+v(426)
       s(5,1)=s(5,1)+v(405)
       s(5,2)=s(5,2)+v(412)
       s(5,3)=s(5,3)+v(418)
       s(5,4)=s(5,4)+v(423)
       s(5,5)=s(5,5)+v(1113)*(v(147)*v(318)+v(149)*v(327))
       s(5,6)=s(5,6)+v(428)
       s(5,7)=s(5,7)+v(429)
       s(5,8)=s(5,8)+v(430)
       s(6,1)=s(6,1)+v(406)
       s(6,2)=s(6,2)+v(413)
       s(6,3)=s(6,3)+v(419)
       s(6,4)=s(6,4)+v(424)
       s(6,5)=s(6,5)+v(428)
       s(6,6)=s(6,6)+v(1113)*(v(147)*v(343)+v(149)*v(360))
       s(6,7)=s(6,7)+v(432)
       s(6,8)=s(6,8)+v(433)
       s(7,1)=s(7,1)+v(407)
       s(7,2)=s(7,2)+v(414)
       s(7,3)=s(7,3)+v(420)
       s(7,4)=s(7,4)+v(425)
       s(7,5)=s(7,5)+v(429)
       s(7,6)=s(7,6)+v(432)
       s(7,7)=s(7,7)+v(1113)*(v(152)*v(320)+v(155)*v(329))
       s(7,8)=s(7,8)+v(435)
       s(8,1)=s(8,1)+v(408)
       s(8,2)=s(8,2)+v(415)
       s(8,3)=s(8,3)+v(421)
       s(8,4)=s(8,4)+v(426)
       s(8,5)=s(8,5)+v(430)
       s(8,6)=s(8,6)+v(433)
       s(8,7)=s(8,7)+v(435)
       s(8,8)=s(8,8)+v(1113)*(v(152)*v(347)+v(155)*v(364))
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
      SUBROUTINE elmt19_ISW05(v,d,xl,ul,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i440
      DOUBLE PRECISION v(1244),d(2),xl(2,4),ul(6,8),m(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i440=1,ngpo
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
      SUBROUTINE elmt19_ISW06(v,d,xl,ul,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i663
      DOUBLE PRECISION v(1244),d(2),xl(2,4),ul(6,8),p(8),ht
     &(0),hp(0),gp(4,4)
      v(731)=ul(1,8)
      v(730)=ul(1,7)
      v(729)=ul(1,6)
      v(728)=ul(1,5)
      v(727)=ul(1,4)
      v(726)=ul(1,3)
      v(725)=ul(1,2)
      v(724)=ul(1,1)
      v(675)=xl(2,4)
      v(674)=xl(1,4)
      v(673)=xl(2,3)
      v(1128)=v(673)-v(675)
      v(672)=xl(1,3)
      v(1130)=v(672)-v(674)
      v(671)=xl(2,2)
      v(1124)=v(671)-v(673)
      v(670)=xl(1,2)
      v(1126)=v(670)-v(672)
      v(669)=xl(2,1)
      v(1127)=v(669)-v(671)
      v(1123)=v(669)-v(675)
      v(668)=xl(1,1)
      v(1129)=v(668)-v(670)
      v(1125)=v(668)-v(674)
      v(665)=d(2)
      v(664)=d(1)
      v(845)=v(664)/(2d0+2d0*v(665))
      v(844)=v(664)/(3d0-6d0*v(665))
      DO i663=1,ngpo
       v(748)=gp(1,i663)
       v(763)=((-1d0)+v(748))/4d0
       v(764)=((-1d0)-v(748))/4d0
       v(769)=v(1123)*v(763)+v(1124)*v(764)
       v(767)=v(1125)*v(763)+v(1126)*v(764)
       v(749)=gp(2,i663)
       v(765)=(1d0+v(749))/4d0
       v(762)=((-1d0)+v(749))/4d0
       v(768)=v(1127)*v(762)+v(1128)*v(765)
       v(766)=v(1129)*v(762)+v(1130)*v(765)
       v(770)=-(v(767)*v(768))+v(766)*v(769)
       v(1122)=gp(4,i663)*v(770)
       v(793)=-(v(769)/v(770))
       v(807)=-(v(765)*v(793))
       v(799)=-(v(762)*v(793))
       v(794)=v(767)/v(770)
       v(810)=-(v(765)*v(794))
       v(801)=-(v(762)*v(794))
       v(795)=-(v(768)/v(770))
       v(808)=v(763)*v(795)
       v(803)=v(764)*v(795)
       v(796)=v(766)/v(770)
       v(811)=v(763)*v(796)
       v(805)=v(764)*v(796)
       v(797)=v(799)+v(808)
       v(798)=v(801)+v(811)
       v(800)=-v(799)+v(803)
       v(802)=-v(801)+v(805)
       v(804)=-v(803)+v(807)
       v(806)=-v(805)+v(810)
       v(809)=-v(807)-v(808)
       v(812)=-v(810)-v(811)
       v(821)=1d0+v(724)*v(797)+v(726)*v(800)+v(728)*v(804)+v(730)*v
     & (809)
       v(822)=v(724)*v(798)+v(726)*v(802)+v(728)*v(806)+v(730)*v(812)
       v(824)=v(725)*v(797)+v(727)*v(800)+v(729)*v(804)+v(731)*v(809)
       v(825)=1d0+v(725)*v(798)+v(727)*v(802)+v(729)*v(806)+v(731)*v
     & (812)
       v(856)=1d0+(v(821)*v(821))+(v(822)*v(822))+(v(824)*v(824))+(v
     & (825)*v(825))
       v(842)=-(v(822)*v(824))+v(821)*v(825)
       v(858)=1d0/v(842)**0.6666666666666666d0
       v(850)=(-3d0)+v(856)*v(858)
       v(854)=(v(845)*(1d0+2d0*d(3)*v(850)+3d0*d(4)*(v(850)*v(850))))
     & /2d0
       v(860)=2d0*v(854)*v(858)
       v(857)=((-1d0)+v(842))*v(844)+((-2d0/3d0)*v(854)*v(856))/v(842
     & )**0.16666666666666669d1
       v(859)=v(825)*v(857)+v(821)*v(860)
       v(861)=-(v(824)*v(857))+v(822)*v(860)
       v(863)=-(v(822)*v(857))+v(824)*v(860)
       v(864)=v(821)*v(857)+v(825)*v(860)
       p(1)=p(1)+v(1122)*(v(797)*v(859)+v(798)*v(861))
       p(2)=p(2)+v(1122)*(v(797)*v(863)+v(798)*v(864))
       p(3)=p(3)+v(1122)*(v(800)*v(859)+v(802)*v(861))
       p(4)=p(4)+v(1122)*(v(800)*v(863)+v(802)*v(864))
       p(5)=p(5)+v(1122)*(v(804)*v(859)+v(806)*v(861))
       p(6)=p(6)+v(1122)*(v(804)*v(863)+v(806)*v(864))
       p(7)=p(7)+v(1122)*(v(809)*v(859)+v(812)*v(861))
       p(8)=p(8)+v(1122)*(v(809)*v(863)+v(812)*v(864))
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt19_ISW09(v,d,xl,ul,c,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i886
      DOUBLE PRECISION v(1244),d(2),xl(2,4),ul(6,8),c(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i886=1,ngpo
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
