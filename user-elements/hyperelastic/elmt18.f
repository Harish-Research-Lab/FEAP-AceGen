!**************************************************************
!* AceGen    7.205 MacOSX (15 Jan 21)                         *
!*           Co. J. Korelc  2020           8 Aug 21 00:33:18  *
!**************************************************************
! User     : Full professional version
! Notebook : neoHookean
! Evaluation time                 : 10 s    Mode  : Optimal
! Number of formulae              : 259     Method: Automatic
! Subroutine                      : elmt18_ISW01 size: 57
! Subroutine                      : elmt18_ISW03 size: 4251
! Subroutine                      : elmt18_ISW05 size: 450
! Subroutine                      : elmt18_ISW06 size: 969
! Subroutine                      : elmt18_ISW09 size: 450
! Total size of Mathematica  code : 6177 subexpressions
! Total size of Fortran code      : 20588 bytes

      subroutine elmt18 (d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

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
      real (kind=8) ::  v(1199)   ! AceGen Storage
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
        utx(1) = "elmt18"

      elseif(isw.eq.1) then                ! Input material set data
        pstyp = ndm                        ! Sets plot dimension (1,2,3)

!       Read number pdes 
        nelu(:)   = 0
        du(:)     = 0
        dofu(:,:) = 0
        call elmt18_ISW01(v,npde,du,nelu,hist1,hist3)

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
          call elmt18_ISW03(v,d,xl,ua,k,c,m,r,hr(nh1),hr(nh2),gp,ngpo)

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
          call elmt18_ISW06(v,d,xl,ua,r,hr(nh1),hr(nh2),gp,ngpo)

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
        call elmt18_ISW05(v,d,xl,ua,m,r,hr(nh1),hr(nh2),gp,ngpo)

      elseif(isw.eq.9) then                ! Compute damping matrix

!       Map dofs
        call SB_ua_set(ul,ndf,nen,ua,du,dofu,nelu,npde)

!       Initialize element arrays
        c(:,:) = 0.0d0
        r(:)   = 0.0d0

!       Call the routine for ISW=9
        call elmt18_ISW09(v,d,xl,ua,c,r,hr(nh1),hr(nh2),gp,ngpo)

      endif

      end subroutine elmt18


!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt18_ISW01(v,npde,du,nelu,hist1,hist3)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER npde,du(10),nelu(10),hist1,hist3
      DOUBLE PRECISION v(1199)
      npde=(1)
      du(1)=(2)
      nelu(1)=(4)
      hist1=(0)
      hist3=(0)
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt18_ISW03(v,d,xl,ul,s,c,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i6
      DOUBLE PRECISION v(1199),d(2),xl(2,4),ul(6,8),s(8,8),c
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
      v(1073)=v(14)-v(16)
      v(13)=xl(1,3)
      v(1075)=v(13)-v(15)
      v(12)=xl(2,2)
      v(1069)=v(12)-v(14)
      v(11)=xl(1,2)
      v(1071)=v(11)-v(13)
      v(10)=xl(2,1)
      v(1072)=v(10)-v(12)
      v(1068)=v(10)-v(16)
      v(9)=xl(1,1)
      v(1074)=-v(11)+v(9)
      v(1070)=-v(15)+v(9)
      v(8)=d(2)
      v(7)=d(1)
      v(186)=v(7)/(2d0+2d0*v(8))
      v(185)=v(7)/(3d0-6d0*v(8))
      DO i6=1,ngpo
       v(89)=gp(1,i6)
       v(104)=((-1d0)+v(89))/4d0
       v(105)=((-1d0)-v(89))/4d0
       v(110)=v(104)*v(1068)+v(105)*v(1069)
       v(108)=v(104)*v(1070)+v(105)*v(1071)
       v(90)=gp(2,i6)
       v(106)=(1d0+v(90))/4d0
       v(103)=((-1d0)+v(90))/4d0
       v(109)=v(103)*v(1072)+v(106)*v(1073)
       v(107)=v(103)*v(1074)+v(106)*v(1075)
       v(111)=-(v(108)*v(109))+v(107)*v(110)
       v(1064)=gp(4,i6)*v(111)
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
       v(1059)=2d0*v(138)
       v(139)=v(142)+v(152)
       v(1060)=2d0*v(139)
       v(141)=-v(140)+v(144)
       v(1057)=2d0*v(141)
       v(143)=-v(142)+v(146)
       v(1058)=2d0*v(143)
       v(145)=-v(144)+v(148)
       v(1055)=2d0*v(145)
       v(147)=-v(146)+v(151)
       v(1056)=2d0*v(147)
       v(150)=-v(148)-v(149)
       v(1052)=2d0*v(150)
       v(153)=-v(151)-v(152)
       v(1053)=2d0*v(153)
       v(162)=1d0+v(138)*v(65)+v(141)*v(67)+v(145)*v(69)+v(150)*v(71)
       v(1067)=v(162)*v(186)
       v(163)=v(139)*v(65)+v(143)*v(67)+v(147)*v(69)+v(153)*v(71)
       v(1066)=v(163)*v(186)
       v(254)=v(153)*v(162)-v(150)*v(163)
       v(252)=v(147)*v(162)-v(145)*v(163)
       v(250)=v(143)*v(162)-v(141)*v(163)
       v(248)=v(139)*v(162)-v(138)*v(163)
       v(165)=v(138)*v(66)+v(141)*v(68)+v(145)*v(70)+v(150)*v(72)
       v(1065)=v(165)*v(186)
       v(166)=1d0+v(139)*v(66)+v(143)*v(68)+v(147)*v(70)+v(153)*v(72)
       v(1063)=v(166)*v(186)
       v(253)=-(v(153)*v(165))+v(150)*v(166)
       v(251)=-(v(147)*v(165))+v(145)*v(166)
       v(249)=-(v(143)*v(165))+v(141)*v(166)
       v(247)=-(v(139)*v(165))+v(138)*v(166)
       v(194)=1d0+(v(162)*v(162))+(v(163)*v(163))+(v(165)*v(165))+(v
     & (166)*v(166))
       v(183)=-(v(163)*v(165))+v(162)*v(166)
       v(1054)=v(185)+((5d0/9d0)*v(186)*v(194))/v(183
     & )**0.26666666666666666d1
       v(255)=1d0/v(183)**0.16666666666666669d1
       v(1061)=(-2d0/3d0)*v(255)
       v(264)=(-1d0/3d0)*(v(186)*v(255))
       v(273)=v(1054)*v(254)+(v(1052)*v(165)+v(1053)*v(166))*v(264)
       v(272)=v(1054)*v(253)+(v(1052)*v(162)+v(1053)*v(163))*v(264)
       v(271)=v(1054)*v(252)+(v(1055)*v(165)+v(1056)*v(166))*v(264)
       v(270)=v(1054)*v(251)+(v(1055)*v(162)+v(1056)*v(163))*v(264)
       v(269)=v(1054)*v(250)+(v(1057)*v(165)+v(1058)*v(166))*v(264)
       v(268)=v(1054)*v(249)+(v(1057)*v(162)+v(1058)*v(163))*v(264)
       v(267)=v(1054)*v(248)+(v(1059)*v(165)+v(1060)*v(166))*v(264)
       v(265)=v(1054)*v(247)+(v(1059)*v(162)+v(1060)*v(163))*v(264)
       v(263)=v(1061)*v(254)
       v(262)=v(1061)*v(253)
       v(261)=v(1061)*v(252)
       v(260)=v(1061)*v(251)
       v(259)=v(1061)*v(250)
       v(258)=v(1061)*v(249)
       v(257)=v(1061)*v(248)
       v(256)=v(1061)*v(247)
       v(1062)=v(186)/v(183)**0.6666666666666666d0
       v(324)=v(1062)*v(153)
       v(325)=v(1063)*v(263)+v(162)*v(273)+v(324)
       v(320)=v(1062)*v(147)
       v(321)=v(1063)*v(261)+v(162)*v(271)+v(320)
       v(316)=v(1062)*v(143)
       v(317)=v(1063)*v(259)+v(162)*v(269)+v(316)
       v(312)=v(1062)*v(139)
       v(307)=v(1062)*v(150)
       v(308)=v(1065)*v(263)-v(163)*v(273)+v(307)
       v(394)=v(1064)*(v(145)*v(308)+v(147)*v(325))
       v(387)=v(1064)*(v(141)*v(308)+v(143)*v(325))
       v(376)=v(1064)*(v(138)*v(308)+v(139)*v(325))
       v(303)=v(1062)*v(145)
       v(304)=v(1065)*v(261)-v(163)*v(271)+v(303)
       v(385)=v(1064)*(v(141)*v(304)+v(143)*v(321))
       v(374)=v(1064)*(v(138)*v(304)+v(139)*v(321))
       v(299)=v(1062)*v(141)
       v(300)=v(1065)*v(259)-v(163)*v(269)+v(299)
       v(372)=v(1064)*(v(138)*v(300)+v(139)*v(317))
       v(295)=v(1062)*v(138)
       v(290)=v(1066)*v(262)-v(165)*v(272)+v(324)
       v(288)=v(1066)*v(260)-v(165)*v(270)+v(320)
       v(286)=v(1066)*v(258)-v(165)*v(268)+v(316)
       v(281)=v(1067)*v(262)+v(166)*v(272)+v(307)
       v(390)=v(1064)*(v(145)*v(281)+v(147)*v(290))
       v(381)=v(1064)*(v(141)*v(281)+v(143)*v(290))
       v(368)=v(1064)*(v(138)*v(281)+v(139)*v(290))
       v(279)=v(1067)*v(260)+v(166)*v(270)+v(303)
       v(379)=v(1064)*(v(141)*v(279)+v(143)*v(288))
       v(366)=v(1064)*(v(138)*v(279)+v(139)*v(288))
       v(277)=v(1067)*v(258)+v(166)*v(268)+v(299)
       v(364)=v(1064)*(v(138)*v(277)+v(139)*v(286))
       v(195)=((-1d0)+v(183))*v(185)+v(194)*v(264)
       v(322)=-(v(150)*v(195))
       v(323)=v(1063)*v(262)+v(162)*v(272)-v(322)
       v(318)=-(v(145)*v(195))
       v(319)=v(1063)*v(260)+v(162)*v(270)-v(318)
       v(314)=-(v(141)*v(195))
       v(305)=v(153)*v(195)
       v(306)=v(1065)*v(262)-v(163)*v(272)-v(305)
       v(393)=v(1064)*(v(145)*v(306)+v(147)*v(323))
       v(386)=v(1064)*(v(141)*v(306)+v(143)*v(323))
       v(375)=v(1064)*(v(138)*v(306)+v(139)*v(323))
       v(301)=v(147)*v(195)
       v(302)=v(1065)*v(260)-v(163)*v(270)-v(301)
       v(384)=v(1064)*(v(141)*v(302)+v(143)*v(319))
       v(373)=v(1064)*(v(138)*v(302)+v(139)*v(319))
       v(297)=v(143)*v(195)
       v(371)=v(1064)*(v(138)*(v(1065)*v(258)-v(163)*v(268)-v(297))+v
     & (139)*(v(1063)*v(258)+v(162)*v(268)-v(314)))
       v(291)=v(1066)*v(263)-v(165)*v(273)+v(322)
       v(289)=v(1066)*v(261)-v(165)*v(271)+v(318)
       v(287)=v(1066)*v(259)-v(165)*v(269)+v(314)
       v(282)=v(1067)*v(263)+v(166)*v(273)+v(305)
       v(396)=v(1064)*(v(150)*v(282)+v(153)*v(291))
       v(391)=v(1064)*(v(145)*v(282)+v(147)*v(291))
       v(382)=v(1064)*(v(141)*v(282)+v(143)*v(291))
       v(369)=v(1064)*(v(138)*v(282)+v(139)*v(291))
       v(280)=v(1067)*v(261)+v(166)*v(271)+v(301)
       v(389)=v(1064)*(v(145)*v(280)+v(147)*v(289))
       v(380)=v(1064)*(v(141)*v(280)+v(143)*v(289))
       v(367)=v(1064)*(v(138)*v(280)+v(139)*v(289))
       v(278)=v(1067)*v(259)+v(166)*v(269)+v(297)
       v(378)=v(1064)*(v(141)*v(278)+v(143)*v(287))
       v(365)=v(1064)*(v(138)*v(278)+v(139)*v(287))
       v(363)=v(1064)*(v(139)*(v(1066)*v(257)-v(165)*v(267))+v(138)*
     & (v(1067)*v(257)+v(166)*v(267)))
       v(197)=v(1062)*v(162)+v(166)*v(195)
       v(198)=v(1062)*v(163)-v(165)*v(195)
       v(200)=v(1062)*v(165)-v(163)*v(195)
       v(201)=v(1062)*v(166)+v(162)*v(195)
       p(1)=p(1)+v(1064)*(v(138)*v(197)+v(139)*v(198))
       p(2)=p(2)+v(1064)*(v(138)*v(200)+v(139)*v(201))
       p(3)=p(3)+v(1064)*(v(141)*v(197)+v(143)*v(198))
       p(4)=p(4)+v(1064)*(v(141)*v(200)+v(143)*v(201))
       p(5)=p(5)+v(1064)*(v(145)*v(197)+v(147)*v(198))
       p(6)=p(6)+v(1064)*(v(145)*v(200)+v(147)*v(201))
       p(7)=p(7)+v(1064)*(v(150)*v(197)+v(153)*v(198))
       p(8)=p(8)+v(1064)*(v(150)*v(200)+v(153)*v(201))
       s(1,1)=s(1,1)+v(1064)*(v(138)*(v(1067)*v(256)+v(166)*v(265)+v
     & (295))+v(139)*(v(1066)*v(256)-v(165)*v(265)+v(312)))
       s(1,2)=s(1,2)+v(363)
       s(1,3)=s(1,3)+v(364)
       s(1,4)=s(1,4)+v(365)
       s(1,5)=s(1,5)+v(366)
       s(1,6)=s(1,6)+v(367)
       s(1,7)=s(1,7)+v(368)
       s(1,8)=s(1,8)+v(369)
       s(2,1)=s(2,1)+v(363)
       s(2,2)=s(2,2)+v(1064)*(v(138)*(v(1065)*v(257)-v(163)*v(267)+v
     & (295))+v(139)*(v(1063)*v(257)+v(162)*v(267)+v(312)))
       s(2,3)=s(2,3)+v(371)
       s(2,4)=s(2,4)+v(372)
       s(2,5)=s(2,5)+v(373)
       s(2,6)=s(2,6)+v(374)
       s(2,7)=s(2,7)+v(375)
       s(2,8)=s(2,8)+v(376)
       s(3,1)=s(3,1)+v(364)
       s(3,2)=s(3,2)+v(371)
       s(3,3)=s(3,3)+v(1064)*(v(141)*v(277)+v(143)*v(286))
       s(3,4)=s(3,4)+v(378)
       s(3,5)=s(3,5)+v(379)
       s(3,6)=s(3,6)+v(380)
       s(3,7)=s(3,7)+v(381)
       s(3,8)=s(3,8)+v(382)
       s(4,1)=s(4,1)+v(365)
       s(4,2)=s(4,2)+v(372)
       s(4,3)=s(4,3)+v(378)
       s(4,4)=s(4,4)+v(1064)*(v(141)*v(300)+v(143)*v(317))
       s(4,5)=s(4,5)+v(384)
       s(4,6)=s(4,6)+v(385)
       s(4,7)=s(4,7)+v(386)
       s(4,8)=s(4,8)+v(387)
       s(5,1)=s(5,1)+v(366)
       s(5,2)=s(5,2)+v(373)
       s(5,3)=s(5,3)+v(379)
       s(5,4)=s(5,4)+v(384)
       s(5,5)=s(5,5)+v(1064)*(v(145)*v(279)+v(147)*v(288))
       s(5,6)=s(5,6)+v(389)
       s(5,7)=s(5,7)+v(390)
       s(5,8)=s(5,8)+v(391)
       s(6,1)=s(6,1)+v(367)
       s(6,2)=s(6,2)+v(374)
       s(6,3)=s(6,3)+v(380)
       s(6,4)=s(6,4)+v(385)
       s(6,5)=s(6,5)+v(389)
       s(6,6)=s(6,6)+v(1064)*(v(145)*v(304)+v(147)*v(321))
       s(6,7)=s(6,7)+v(393)
       s(6,8)=s(6,8)+v(394)
       s(7,1)=s(7,1)+v(368)
       s(7,2)=s(7,2)+v(375)
       s(7,3)=s(7,3)+v(381)
       s(7,4)=s(7,4)+v(386)
       s(7,5)=s(7,5)+v(390)
       s(7,6)=s(7,6)+v(393)
       s(7,7)=s(7,7)+v(1064)*(v(150)*v(281)+v(153)*v(290))
       s(7,8)=s(7,8)+v(396)
       s(8,1)=s(8,1)+v(369)
       s(8,2)=s(8,2)+v(376)
       s(8,3)=s(8,3)+v(382)
       s(8,4)=s(8,4)+v(387)
       s(8,5)=s(8,5)+v(391)
       s(8,6)=s(8,6)+v(394)
       s(8,7)=s(8,7)+v(396)
       s(8,8)=s(8,8)+v(1064)*(v(150)*v(308)+v(153)*v(325))
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
      SUBROUTINE elmt18_ISW05(v,d,xl,ul,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i401
      DOUBLE PRECISION v(1199),d(2),xl(2,4),ul(6,8),m(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i401=1,ngpo
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
      SUBROUTINE elmt18_ISW06(v,d,xl,ul,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i618
      DOUBLE PRECISION v(1199),d(2),xl(2,4),ul(6,8),p(8),ht
     &(0),hp(0),gp(4,4)
      v(684)=ul(1,8)
      v(683)=ul(1,7)
      v(682)=ul(1,6)
      v(681)=ul(1,5)
      v(680)=ul(1,4)
      v(679)=ul(1,3)
      v(678)=ul(1,2)
      v(677)=ul(1,1)
      v(628)=xl(2,4)
      v(627)=xl(1,4)
      v(626)=xl(2,3)
      v(1083)=v(626)-v(628)
      v(625)=xl(1,3)
      v(1085)=v(625)-v(627)
      v(624)=xl(2,2)
      v(1079)=v(624)-v(626)
      v(623)=xl(1,2)
      v(1081)=v(623)-v(625)
      v(622)=xl(2,1)
      v(1082)=v(622)-v(624)
      v(1078)=v(622)-v(628)
      v(621)=xl(1,1)
      v(1084)=v(621)-v(623)
      v(1080)=v(621)-v(627)
      v(620)=d(2)
      v(619)=d(1)
      v(798)=v(619)/(2d0+2d0*v(620))
      v(797)=v(619)/(3d0-6d0*v(620))
      DO i618=1,ngpo
       v(701)=gp(1,i618)
       v(716)=((-1d0)+v(701))/4d0
       v(717)=((-1d0)-v(701))/4d0
       v(722)=v(1078)*v(716)+v(1079)*v(717)
       v(720)=v(1080)*v(716)+v(1081)*v(717)
       v(702)=gp(2,i618)
       v(718)=(1d0+v(702))/4d0
       v(715)=((-1d0)+v(702))/4d0
       v(721)=v(1082)*v(715)+v(1083)*v(718)
       v(719)=v(1084)*v(715)+v(1085)*v(718)
       v(723)=-(v(720)*v(721))+v(719)*v(722)
       v(1077)=gp(4,i618)*v(723)
       v(746)=-(v(722)/v(723))
       v(760)=-(v(718)*v(746))
       v(752)=-(v(715)*v(746))
       v(747)=v(720)/v(723)
       v(763)=-(v(718)*v(747))
       v(754)=-(v(715)*v(747))
       v(748)=-(v(721)/v(723))
       v(761)=v(716)*v(748)
       v(756)=v(717)*v(748)
       v(749)=v(719)/v(723)
       v(764)=v(716)*v(749)
       v(758)=v(717)*v(749)
       v(750)=v(752)+v(761)
       v(751)=v(754)+v(764)
       v(753)=-v(752)+v(756)
       v(755)=-v(754)+v(758)
       v(757)=-v(756)+v(760)
       v(759)=-v(758)+v(763)
       v(762)=-v(760)-v(761)
       v(765)=-v(763)-v(764)
       v(774)=1d0+v(677)*v(750)+v(679)*v(753)+v(681)*v(757)+v(683)*v
     & (762)
       v(775)=v(677)*v(751)+v(679)*v(755)+v(681)*v(759)+v(683)*v(765)
       v(777)=v(678)*v(750)+v(680)*v(753)+v(682)*v(757)+v(684)*v(762)
       v(778)=1d0+v(678)*v(751)+v(680)*v(755)+v(682)*v(759)+v(684)*v
     & (765)
       v(795)=-(v(775)*v(777))+v(774)*v(778)
       v(1076)=v(798)/v(795)**0.6666666666666666d0
       v(807)=((-1d0)+v(795))*v(797)-((1d0+(v(774)*v(774))+(v(775)*v
     & (775))+(v(777)*v(777))+(v(778)*v(778)))*v(798))/(3d0*v(795
     & )**0.16666666666666669d1)
       v(809)=v(1076)*v(774)+v(778)*v(807)
       v(810)=v(1076)*v(775)-v(777)*v(807)
       v(812)=v(1076)*v(777)-v(775)*v(807)
       v(813)=v(1076)*v(778)+v(774)*v(807)
       p(1)=p(1)+v(1077)*(v(750)*v(809)+v(751)*v(810))
       p(2)=p(2)+v(1077)*(v(750)*v(812)+v(751)*v(813))
       p(3)=p(3)+v(1077)*(v(753)*v(809)+v(755)*v(810))
       p(4)=p(4)+v(1077)*(v(753)*v(812)+v(755)*v(813))
       p(5)=p(5)+v(1077)*(v(757)*v(809)+v(759)*v(810))
       p(6)=p(6)+v(1077)*(v(757)*v(812)+v(759)*v(813))
       p(7)=p(7)+v(1077)*(v(762)*v(809)+v(765)*v(810))
       p(8)=p(8)+v(1077)*(v(762)*v(812)+v(765)*v(813))
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt18_ISW09(v,d,xl,ul,c,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i835
      DOUBLE PRECISION v(1199),d(2),xl(2,4),ul(6,8),c(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i835=1,ngpo
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
