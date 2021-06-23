!**************************************************************
!* AceGen    7.205 MacOSX (15 Jan 21)                         *
!*           Co. J. Korelc  2020           21 Jun 21 15:46:00 *
!**************************************************************
! User     : Full professional version
! Notebook : neoHookean
! Evaluation time                 : 8 s     Mode  : Optimal
! Number of formulae              : 243     Method: Automatic
! Subroutine                      : elmt18_ISW01 size: 57
! Subroutine                      : elmt18_ISW03 size: 3999
! Subroutine                      : elmt18_ISW05 size: 450
! Subroutine                      : elmt18_ISW06 size: 968
! Subroutine                      : elmt18_ISW09 size: 450
! Total size of Mathematica  code : 5924 subexpressions
! Total size of Fortran code      : 19524 bytes

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
      real (kind=8) ::  v(1088)   ! AceGen Storage
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

      elseif(isw.eq.9) then                ! Compute mass matrix

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
      DOUBLE PRECISION v(1088)
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
      DOUBLE PRECISION v(1088),d(2),xl(2,4),ul(6,8),s(8,8),c
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
      v(961)=v(14)-v(16)
      v(13)=xl(1,3)
      v(963)=v(13)-v(15)
      v(12)=xl(2,2)
      v(957)=v(12)-v(14)
      v(11)=xl(1,2)
      v(959)=v(11)-v(13)
      v(10)=xl(2,1)
      v(960)=v(10)-v(12)
      v(956)=v(10)-v(16)
      v(9)=xl(1,1)
      v(962)=-v(11)+v(9)
      v(958)=-v(15)+v(9)
      v(8)=d(2)
      v(934)=d(1)/(1d0+v(8))
      v(177)=(v(8)*v(934))/(1d0-2d0*v(8))
      v(175)=v(934)/2d0
      DO i6=1,ngpo
       v(89)=gp(1,i6)
       v(104)=((-1d0)+v(89))/4d0
       v(105)=((-1d0)-v(89))/4d0
       v(110)=v(104)*v(956)+v(105)*v(957)
       v(108)=v(104)*v(958)+v(105)*v(959)
       v(90)=gp(2,i6)
       v(106)=(1d0+v(90))/4d0
       v(103)=((-1d0)+v(90))/4d0
       v(109)=v(103)*v(960)+v(106)*v(961)
       v(107)=v(103)*v(962)+v(106)*v(963)
       v(111)=-(v(108)*v(109))+v(107)*v(110)
       v(938)=gp(4,i6)*v(111)
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
       v(139)=v(142)+v(152)
       v(141)=-v(140)+v(144)
       v(944)=v(141)*v(175)
       v(143)=-v(142)+v(146)
       v(948)=v(143)*v(175)
       v(145)=-v(144)+v(148)
       v(942)=v(145)*v(175)
       v(147)=-v(146)+v(151)
       v(946)=v(147)*v(175)
       v(150)=-v(148)-v(149)
       v(941)=v(150)*v(175)
       v(153)=-v(151)-v(152)
       v(945)=v(153)*v(175)
       v(158)=1d0+v(138)*v(65)+v(141)*v(67)+v(145)*v(69)+v(150)*v(71)
       v(954)=v(158)*v(177)
       v(159)=v(139)*v(65)+v(143)*v(67)+v(147)*v(69)+v(153)*v(71)
       v(953)=-(v(159)*v(177))
       v(212)=v(153)*v(158)-v(150)*v(159)
       v(210)=v(147)*v(158)-v(145)*v(159)
       v(208)=v(143)*v(158)-v(141)*v(159)
       v(206)=v(139)*v(158)-v(138)*v(159)
       v(160)=v(138)*v(66)+v(141)*v(68)+v(145)*v(70)+v(150)*v(72)
       v(951)=-(v(160)*v(177))
       v(161)=1d0+v(139)*v(66)+v(143)*v(68)+v(147)*v(70)+v(153)*v(72)
       v(950)=v(161)*v(177)
       v(211)=-(v(153)*v(160))+v(150)*v(161)
       v(209)=-(v(147)*v(160))+v(145)*v(161)
       v(207)=-(v(143)*v(160))+v(141)*v(161)
       v(205)=-(v(139)*v(160))+v(138)*v(161)
       v(952)=-(v(175)*v(205))
       v(167)=-(v(159)*v(160))+v(158)*v(161)
       v(214)=1d0/v(167)**2
       v(955)=v(205)*v(214)
       v(939)=v(158)*v(214)
       v(937)=v(177)+v(175)*v(214)
       v(947)=v(158)*v(937)
       v(943)=-(v(159)*v(937))
       v(940)=-(v(160)*v(937))
       v(936)=v(161)*v(937)
       v(935)=-(v(159)*v(214))
       v(260)=v(211)*v(935)
       v(258)=v(211)*v(936)+v(941)
       v(256)=v(209)*v(935)
       v(254)=v(209)*v(936)+v(942)
       v(252)=v(207)*v(935)
       v(250)=v(207)*v(936)+v(944)
       v(244)=v(211)*v(939)
       v(242)=v(211)*v(940)+v(945)
       v(326)=(v(147)*v(242)+v(145)*v(258))*v(938)
       v(317)=(v(143)*v(242)+v(141)*v(258))*v(938)
       v(304)=(v(139)*v(242)+v(138)*v(258))*v(938)
       v(240)=v(209)*v(939)
       v(238)=v(209)*v(940)+v(946)
       v(315)=(v(143)*v(238)+v(141)*v(254))*v(938)
       v(302)=(v(139)*v(238)+v(138)*v(254))*v(938)
       v(236)=v(207)*v(939)
       v(234)=v(207)*v(940)+v(948)
       v(300)=(v(139)*v(234)+v(138)*v(250))*v(938)
       v(229)=v(941)+v(212)*v(943)
       v(227)=v(942)+v(210)*v(943)
       v(225)=v(208)*v(943)+v(944)
       v(221)=v(945)+v(212)*v(947)
       v(330)=(v(147)*v(221)+v(145)*v(229))*v(938)
       v(323)=(v(143)*v(221)+v(141)*v(229))*v(938)
       v(312)=(v(139)*v(221)+v(138)*v(229))*v(938)
       v(219)=v(946)+v(210)*v(947)
       v(321)=(v(143)*v(219)+v(141)*v(227))*v(938)
       v(310)=(v(139)*v(219)+v(138)*v(227))*v(938)
       v(217)=v(208)*v(947)+v(948)
       v(308)=(v(139)*v(217)+v(138)*v(225))*v(938)
       v(183)=(-1d0)+v(167)
       v(949)=-(v(177)*v(183))
       v(259)=v(153)*v(949)
       v(261)=-v(259)+v(175)*v(260)+v(212)*v(950)
       v(255)=v(147)*v(949)
       v(257)=-v(255)+v(175)*v(256)+v(210)*v(950)
       v(251)=v(143)*v(949)
       v(253)=-v(251)+v(175)*v(252)+v(208)*v(950)
       v(243)=-(v(150)*v(949))
       v(245)=-v(243)+v(175)*v(244)+v(212)*v(951)
       v(332)=(v(153)*v(245)+v(150)*v(261))*v(938)
       v(327)=(v(147)*v(245)+v(145)*v(261))*v(938)
       v(318)=(v(143)*v(245)+v(141)*v(261))*v(938)
       v(305)=(v(139)*v(245)+v(138)*v(261))*v(938)
       v(239)=-(v(145)*v(949))
       v(241)=-v(239)+v(175)*v(240)+v(210)*v(951)
       v(325)=(v(147)*v(241)+v(145)*v(257))*v(938)
       v(316)=(v(143)*v(241)+v(141)*v(257))*v(938)
       v(303)=(v(139)*v(241)+v(138)*v(257))*v(938)
       v(235)=-(v(141)*v(949))
       v(237)=-v(235)+v(175)*v(236)+v(208)*v(951)
       v(314)=(v(143)*v(237)+v(141)*v(253))*v(938)
       v(301)=(v(139)*v(237)+v(138)*v(253))*v(938)
       v(299)=v(938)*(-(v(138)*(-(v(206)*v(950))+v(935)*v(952)))-v
     & (139)*(-(v(206)*v(951))+v(939)*v(952)))
       v(228)=v(259)+v(175)*(v(153)/v(167)+v(260))+v(211)*v(953)
       v(226)=v(255)+v(175)*(v(147)/v(167)+v(256))+v(209)*v(953)
       v(220)=v(243)+v(175)*(-(v(150)/v(167))+v(244))+v(211)*v(954)
       v(329)=(v(147)*v(220)+v(145)*v(228))*v(938)
       v(322)=(v(143)*v(220)+v(141)*v(228))*v(938)
       v(311)=(v(139)*v(220)+v(138)*v(228))*v(938)
       v(218)=v(239)+v(175)*(-(v(145)/v(167))+v(240))+v(209)*v(954)
       v(320)=(v(143)*v(218)+v(141)*v(226))*v(938)
       v(309)=(v(139)*v(218)+v(138)*v(226))*v(938)
       v(307)=v(938)*(v(138)*(v(251)+v(175)*(v(143)/v(167)+v(252))+v
     & (207)*v(953))+v(139)*(v(235)+v(175)*(-(v(141)/v(167))+v(236))
     & +v(207)*v(954)))
       v(187)=(v(161)-v(158)/v(167))*v(175)-v(158)*v(949)
       v(186)=(v(160)+v(159)/v(167))*v(175)+v(183)*v(953)
       v(185)=(v(159)+v(160)/v(167))*v(175)+v(183)*v(951)
       v(184)=(v(158)-v(161)/v(167))*v(175)-v(161)*v(949)
       p(1)=p(1)+(v(138)*v(184)+v(139)*v(185))*v(938)
       p(2)=p(2)+(v(138)*v(186)+v(139)*v(187))*v(938)
       p(3)=p(3)+(v(141)*v(184)+v(143)*v(185))*v(938)
       p(4)=p(4)+(v(141)*v(186)+v(143)*v(187))*v(938)
       p(5)=p(5)+(v(145)*v(184)+v(147)*v(185))*v(938)
       p(6)=p(6)+(v(145)*v(186)+v(147)*v(187))*v(938)
       p(7)=p(7)+(v(150)*v(184)+v(153)*v(185))*v(938)
       p(8)=p(8)+(v(150)*v(186)+v(153)*v(187))*v(938)
       s(1,1)=s(1,1)+v(938)*(v(139)*(v(205)*v(951)+v(175)*(v(139)-v
     & (160)*v(955)))+v(138)*(v(205)*v(950)+v(175)*(v(138)+v(161)*v
     & (955))))
       s(1,2)=s(1,2)+v(299)
       s(1,3)=s(1,3)+v(300)
       s(1,4)=s(1,4)+v(301)
       s(1,5)=s(1,5)+v(302)
       s(1,6)=s(1,6)+v(303)
       s(1,7)=s(1,7)+v(304)
       s(1,8)=s(1,8)+v(305)
       s(2,1)=s(2,1)+v(299)
       s(2,2)=s(2,2)+v(938)*(v(138)*(v(175)*(v(138)+v(206)*v(935))+v
     & (206)*v(953))+v(139)*(v(175)*(v(139)+v(206)*v(939))+v(206)*v
     & (954)))
       s(2,3)=s(2,3)+v(307)
       s(2,4)=s(2,4)+v(308)
       s(2,5)=s(2,5)+v(309)
       s(2,6)=s(2,6)+v(310)
       s(2,7)=s(2,7)+v(311)
       s(2,8)=s(2,8)+v(312)
       s(3,1)=s(3,1)+v(300)
       s(3,2)=s(3,2)+v(307)
       s(3,3)=s(3,3)+(v(143)*v(234)+v(141)*v(250))*v(938)
       s(3,4)=s(3,4)+v(314)
       s(3,5)=s(3,5)+v(315)
       s(3,6)=s(3,6)+v(316)
       s(3,7)=s(3,7)+v(317)
       s(3,8)=s(3,8)+v(318)
       s(4,1)=s(4,1)+v(301)
       s(4,2)=s(4,2)+v(308)
       s(4,3)=s(4,3)+v(314)
       s(4,4)=s(4,4)+(v(143)*v(217)+v(141)*v(225))*v(938)
       s(4,5)=s(4,5)+v(320)
       s(4,6)=s(4,6)+v(321)
       s(4,7)=s(4,7)+v(322)
       s(4,8)=s(4,8)+v(323)
       s(5,1)=s(5,1)+v(302)
       s(5,2)=s(5,2)+v(309)
       s(5,3)=s(5,3)+v(315)
       s(5,4)=s(5,4)+v(320)
       s(5,5)=s(5,5)+(v(147)*v(238)+v(145)*v(254))*v(938)
       s(5,6)=s(5,6)+v(325)
       s(5,7)=s(5,7)+v(326)
       s(5,8)=s(5,8)+v(327)
       s(6,1)=s(6,1)+v(303)
       s(6,2)=s(6,2)+v(310)
       s(6,3)=s(6,3)+v(316)
       s(6,4)=s(6,4)+v(321)
       s(6,5)=s(6,5)+v(325)
       s(6,6)=s(6,6)+(v(147)*v(219)+v(145)*v(227))*v(938)
       s(6,7)=s(6,7)+v(329)
       s(6,8)=s(6,8)+v(330)
       s(7,1)=s(7,1)+v(304)
       s(7,2)=s(7,2)+v(311)
       s(7,3)=s(7,3)+v(317)
       s(7,4)=s(7,4)+v(322)
       s(7,5)=s(7,5)+v(326)
       s(7,6)=s(7,6)+v(329)
       s(7,7)=s(7,7)+(v(153)*v(242)+v(150)*v(258))*v(938)
       s(7,8)=s(7,8)+v(332)
       s(8,1)=s(8,1)+v(305)
       s(8,2)=s(8,2)+v(312)
       s(8,3)=s(8,3)+v(318)
       s(8,4)=s(8,4)+v(323)
       s(8,5)=s(8,5)+v(327)
       s(8,6)=s(8,6)+v(330)
       s(8,7)=s(8,7)+v(332)
       s(8,8)=s(8,8)+(v(153)*v(221)+v(150)*v(229))*v(938)
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
      INTEGER ngpo,i337
      DOUBLE PRECISION v(1088),d(2),xl(2,4),ul(6,8),m(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i337=1,ngpo
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
      INTEGER ngpo,i536
      DOUBLE PRECISION v(1088),d(2),xl(2,4),ul(6,8),p(8),ht
     &(0),hp(0),gp(4,4)
      v(602)=ul(1,8)
      v(601)=ul(1,7)
      v(600)=ul(1,6)
      v(599)=ul(1,5)
      v(598)=ul(1,4)
      v(597)=ul(1,3)
      v(596)=ul(1,2)
      v(595)=ul(1,1)
      v(546)=xl(2,4)
      v(545)=xl(1,4)
      v(544)=xl(2,3)
      v(972)=v(544)-v(546)
      v(543)=xl(1,3)
      v(974)=v(543)-v(545)
      v(542)=xl(2,2)
      v(968)=v(542)-v(544)
      v(541)=xl(1,2)
      v(970)=v(541)-v(543)
      v(540)=xl(2,1)
      v(971)=v(540)-v(542)
      v(967)=v(540)-v(546)
      v(539)=xl(1,1)
      v(973)=v(539)-v(541)
      v(969)=v(539)-v(545)
      v(538)=d(2)
      v(964)=d(1)/(1d0+v(538))
      v(707)=(v(538)*v(964))/(1d0-2d0*v(538))
      v(705)=v(964)/2d0
      DO i536=1,ngpo
       v(619)=gp(1,i536)
       v(634)=((-1d0)+v(619))/4d0
       v(635)=((-1d0)-v(619))/4d0
       v(640)=v(634)*v(967)+v(635)*v(968)
       v(638)=v(634)*v(969)+v(635)*v(970)
       v(620)=gp(2,i536)
       v(636)=(1d0+v(620))/4d0
       v(633)=((-1d0)+v(620))/4d0
       v(639)=v(633)*v(971)+v(636)*v(972)
       v(637)=v(633)*v(973)+v(636)*v(974)
       v(641)=-(v(638)*v(639))+v(637)*v(640)
       v(966)=gp(4,i536)*v(641)
       v(664)=-(v(640)/v(641))
       v(678)=-(v(636)*v(664))
       v(670)=-(v(633)*v(664))
       v(665)=v(638)/v(641)
       v(681)=-(v(636)*v(665))
       v(672)=-(v(633)*v(665))
       v(666)=-(v(639)/v(641))
       v(679)=v(634)*v(666)
       v(674)=v(635)*v(666)
       v(667)=v(637)/v(641)
       v(682)=v(634)*v(667)
       v(676)=v(635)*v(667)
       v(668)=v(670)+v(679)
       v(669)=v(672)+v(682)
       v(671)=-v(670)+v(674)
       v(673)=-v(672)+v(676)
       v(675)=-v(674)+v(678)
       v(677)=-v(676)+v(681)
       v(680)=-v(678)-v(679)
       v(683)=-v(681)-v(682)
       v(688)=1d0+v(595)*v(668)+v(597)*v(671)+v(599)*v(675)+v(601)*v
     & (680)
       v(689)=v(595)*v(669)+v(597)*v(673)+v(599)*v(677)+v(601)*v(683)
       v(690)=v(596)*v(668)+v(598)*v(671)+v(600)*v(675)+v(602)*v(680)
       v(691)=1d0+v(596)*v(669)+v(598)*v(673)+v(600)*v(677)+v(602)*v
     & (683)
       v(697)=-(v(689)*v(690))+v(688)*v(691)
       v(965)=((-1d0)+v(697))*v(707)
       v(717)=(v(691)-v(688)/v(697))*v(705)+v(688)*v(965)
       v(716)=(v(690)+v(689)/v(697))*v(705)-v(689)*v(965)
       v(715)=(v(689)+v(690)/v(697))*v(705)-v(690)*v(965)
       v(714)=(v(688)-v(691)/v(697))*v(705)+v(691)*v(965)
       p(1)=p(1)+(v(668)*v(714)+v(669)*v(715))*v(966)
       p(2)=p(2)+(v(668)*v(716)+v(669)*v(717))*v(966)
       p(3)=p(3)+(v(671)*v(714)+v(673)*v(715))*v(966)
       p(4)=p(4)+(v(671)*v(716)+v(673)*v(717))*v(966)
       p(5)=p(5)+(v(675)*v(714)+v(677)*v(715))*v(966)
       p(6)=p(6)+(v(675)*v(716)+v(677)*v(717))*v(966)
       p(7)=p(7)+(v(680)*v(714)+v(683)*v(715))*v(966)
       p(8)=p(8)+(v(680)*v(716)+v(683)*v(717))*v(966)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt18_ISW09(v,d,xl,ul,c,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i735
      DOUBLE PRECISION v(1088),d(2),xl(2,4),ul(6,8),c(8,8),p
     &(8),ht(0),hp(0),gp(4,4)
      DO i735=1,ngpo
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
