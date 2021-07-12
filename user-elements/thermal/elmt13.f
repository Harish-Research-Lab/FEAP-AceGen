!**************************************************************
!* AceGen    7.205 MacOSX (15 Jan 21)                         *
!*           Co. J. Korelc  2020           12 Jul 21 01:32:25 *
!**************************************************************
! User     : Full professional version
! Notebook : thermalSteady
! Evaluation time                 : 4 s     Mode  : Optimal
! Number of formulae              : 126     Method: Automatic
! Subroutine                      : elmt13_ISW01 size: 57
! Subroutine                      : elmt13_ISW03 size: 1159
! Subroutine                      : elmt13_ISW05 size: 114
! Subroutine                      : elmt13_ISW06 size: 623
! Subroutine                      : elmt13_ISW09 size: 114
! Total size of Mathematica  code : 2067 subexpressions
! Total size of Fortran code      : 7627 bytes

      subroutine elmt13 (d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

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
      real (kind=8) ::  v(592)   ! AceGen Storage
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
        utx(1) = "elmt13"

      elseif(isw.eq.1) then                ! Input material set data
        pstyp = ndm                        ! Sets plot dimension (1,2,3)

!       Read number pdes 
        nelu(:)   = 0
        du(:)     = 0
        dofu(:,:) = 0
        call elmt13_ISW01(v,npde,du,nelu,hist1,hist3)

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
        
        datades(1)="rho0 -Density"

        datades(2)="kappa -Thermal conductivity"

        datades(3)="Cv -Heat capacity"

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
          call elmt13_ISW03(v,d,xl,ua,k,c,m,r,hr(nh1),hr(nh2),gp,ngpo)

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
          call elmt13_ISW06(v,d,xl,ua,r,hr(nh1),hr(nh2),gp,ngpo)

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
        call elmt13_ISW05(v,d,xl,ua,m,r,hr(nh1),hr(nh2),gp,ngpo)

      elseif(isw.eq.9) then                ! Compute damping matrix

!       Map dofs
        call SB_ua_set(ul,ndf,nen,ua,du,dofu,nelu,npde)

!       Initialize element arrays
        c(:,:) = 0.0d0
        r(:)   = 0.0d0

!       Call the routine for ISW=9
        call elmt13_ISW09(v,d,xl,ua,c,r,hr(nh1),hr(nh2),gp,ngpo)

      endif

      end subroutine elmt13


!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt13_ISW01(v,npde,du,nelu,hist1,hist3)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER npde,du(10),nelu(10),hist1,hist3
      DOUBLE PRECISION v(592)
      npde=(1)
      du(1)=(1)
      nelu(1)=(4)
      hist1=(0)
      hist3=(0)
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt13_ISW03(v,d,xl,ul,s,c,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i6
      DOUBLE PRECISION v(592),d(3),xl(2,4),ul(6,4),s(4,4),c
     &(4,4),m(4,4),p(4),ht(0),hp(0),gp(4,4)
      v(45)=ul(1,4)
      v(44)=ul(1,3)
      v(43)=ul(1,2)
      v(42)=ul(1,1)
      v(17)=xl(2,4)
      v(16)=xl(1,4)
      v(15)=xl(2,3)
      v(467)=v(15)-v(17)
      v(14)=xl(1,3)
      v(469)=v(14)-v(16)
      v(13)=xl(2,2)
      v(463)=v(13)-v(15)
      v(12)=xl(1,2)
      v(465)=v(12)-v(14)
      v(11)=xl(2,1)
      v(466)=v(11)-v(13)
      v(462)=v(11)-v(17)
      v(10)=xl(1,1)
      v(468)=v(10)-v(12)
      v(464)=v(10)-v(16)
      DO i6=1,ngpo
       v(54)=gp(1,i6)
       v(69)=((-1d0)+v(54))/4d0
       v(70)=((-1d0)-v(54))/4d0
       v(75)=v(462)*v(69)+v(463)*v(70)
       v(73)=v(464)*v(69)+v(465)*v(70)
       v(55)=gp(2,i6)
       v(71)=(1d0+v(55))/4d0
       v(68)=((-1d0)+v(55))/4d0
       v(74)=v(466)*v(68)+v(467)*v(71)
       v(72)=v(468)*v(68)+v(469)*v(71)
       v(76)=-(v(73)*v(74))+v(72)*v(75)
       v(461)=d(2)*gp(4,i6)*v(76)
       v(82)=-(v(75)/v(76))
       v(96)=-(v(71)*v(82))
       v(88)=-(v(68)*v(82))
       v(83)=v(73)/v(76)
       v(99)=-(v(71)*v(83))
       v(90)=-(v(68)*v(83))
       v(84)=-(v(74)/v(76))
       v(97)=v(69)*v(84)
       v(92)=v(70)*v(84)
       v(85)=v(72)/v(76)
       v(100)=v(69)*v(85)
       v(94)=v(70)*v(85)
       v(86)=v(88)+v(97)
       v(87)=v(100)+v(90)
       v(89)=-v(88)+v(92)
       v(91)=-v(90)+v(94)
       v(125)=v(461)*(v(86)*v(89)+v(87)*v(91))
       v(93)=-v(92)+v(96)
       v(95)=-v(94)+v(99)
       v(129)=v(461)*(v(89)*v(93)+v(91)*v(95))
       v(126)=v(461)*(v(86)*v(93)+v(87)*v(95))
       v(98)=-v(96)-v(97)
       v(102)=v(42)*v(86)+v(43)*v(89)+v(44)*v(93)+v(45)*v(98)
       v(101)=-v(100)-v(99)
       v(132)=v(461)*(v(101)*v(95)+v(93)*v(98))
       v(130)=v(461)*(v(101)*v(91)+v(89)*v(98))
       v(127)=v(461)*(v(101)*v(87)+v(86)*v(98))
       v(103)=v(101)*v(45)+v(42)*v(87)+v(43)*v(91)+v(44)*v(95)
       p(1)=p(1)+v(461)*(v(102)*v(86)+v(103)*v(87))
       p(2)=p(2)+v(461)*(v(102)*v(89)+v(103)*v(91))
       p(3)=p(3)+v(461)*(v(102)*v(93)+v(103)*v(95))
       p(4)=p(4)+v(461)*(v(101)*v(103)+v(102)*v(98))
       s(1,1)=s(1,1)+v(461)*((v(86)*v(86))+(v(87)*v(87)))
       s(1,2)=s(1,2)+v(125)
       s(1,3)=s(1,3)+v(126)
       s(1,4)=s(1,4)+v(127)
       s(2,1)=s(2,1)+v(125)
       s(2,2)=s(2,2)+v(461)*((v(89)*v(89))+(v(91)*v(91)))
       s(2,3)=s(2,3)+v(129)
       s(2,4)=s(2,4)+v(130)
       s(3,1)=s(3,1)+v(126)
       s(3,2)=s(3,2)+v(129)
       s(3,3)=s(3,3)+v(461)*((v(93)*v(93))+(v(95)*v(95)))
       s(3,4)=s(3,4)+v(132)
       s(4,1)=s(4,1)+v(127)
       s(4,2)=s(4,2)+v(130)
       s(4,3)=s(4,3)+v(132)
       s(4,4)=s(4,4)+v(461)*((v(101)*v(101))+(v(98)*v(98)))
       c(1,1)=0d0+c(1,1)
       c(1,2)=0d0+c(1,2)
       c(1,3)=0d0+c(1,3)
       c(1,4)=0d0+c(1,4)
       c(2,1)=0d0+c(2,1)
       c(2,2)=0d0+c(2,2)
       c(2,3)=0d0+c(2,3)
       c(2,4)=0d0+c(2,4)
       c(3,1)=0d0+c(3,1)
       c(3,2)=0d0+c(3,2)
       c(3,3)=0d0+c(3,3)
       c(3,4)=0d0+c(3,4)
       c(4,1)=0d0+c(4,1)
       c(4,2)=0d0+c(4,2)
       c(4,3)=0d0+c(4,3)
       c(4,4)=0d0+c(4,4)
       m(1,1)=0d0+m(1,1)
       m(1,2)=0d0+m(1,2)
       m(1,3)=0d0+m(1,3)
       m(1,4)=0d0+m(1,4)
       m(2,1)=0d0+m(2,1)
       m(2,2)=0d0+m(2,2)
       m(2,3)=0d0+m(2,3)
       m(2,4)=0d0+m(2,4)
       m(3,1)=0d0+m(3,1)
       m(3,2)=0d0+m(3,2)
       m(3,3)=0d0+m(3,3)
       m(3,4)=0d0+m(3,4)
       m(4,1)=0d0+m(4,1)
       m(4,2)=0d0+m(4,2)
       m(4,3)=0d0+m(4,3)
       m(4,4)=0d0+m(4,4)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt13_ISW05(v,d,xl,ul,m,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i137
      DOUBLE PRECISION v(592),d(3),xl(2,4),ul(6,4),m(4,4),p
     &(4),ht(0),hp(0),gp(4,4)
      DO i137=1,ngpo
       m(1,1)=0d0+m(1,1)
       m(1,2)=0d0+m(1,2)
       m(1,3)=0d0+m(1,3)
       m(1,4)=0d0+m(1,4)
       m(2,1)=0d0+m(2,1)
       m(2,2)=0d0+m(2,2)
       m(2,3)=0d0+m(2,3)
       m(2,4)=0d0+m(2,4)
       m(3,1)=0d0+m(3,1)
       m(3,2)=0d0+m(3,2)
       m(3,3)=0d0+m(3,3)
       m(3,4)=0d0+m(3,4)
       m(4,1)=0d0+m(4,1)
       m(4,2)=0d0+m(4,2)
       m(4,3)=0d0+m(4,3)
       m(4,4)=0d0+m(4,4)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt13_ISW06(v,d,xl,ul,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i245
      DOUBLE PRECISION v(592),d(3),xl(2,4),ul(6,4),p(4),ht
     &(0),hp(0),gp(4,4)
      v(284)=ul(1,4)
      v(283)=ul(1,3)
      v(282)=ul(1,2)
      v(281)=ul(1,1)
      v(256)=xl(2,4)
      v(255)=xl(1,4)
      v(254)=xl(2,3)
      v(476)=v(254)-v(256)
      v(253)=xl(1,3)
      v(478)=v(253)-v(255)
      v(252)=xl(2,2)
      v(472)=v(252)-v(254)
      v(251)=xl(1,2)
      v(474)=v(251)-v(253)
      v(250)=xl(2,1)
      v(475)=v(250)-v(252)
      v(471)=v(250)-v(256)
      v(249)=xl(1,1)
      v(477)=v(249)-v(251)
      v(473)=v(249)-v(255)
      DO i245=1,ngpo
       v(293)=gp(1,i245)
       v(308)=((-1d0)+v(293))/4d0
       v(309)=((-1d0)-v(293))/4d0
       v(314)=v(308)*v(471)+v(309)*v(472)
       v(312)=v(308)*v(473)+v(309)*v(474)
       v(294)=gp(2,i245)
       v(310)=(1d0+v(294))/4d0
       v(307)=((-1d0)+v(294))/4d0
       v(313)=v(307)*v(475)+v(310)*v(476)
       v(311)=v(307)*v(477)+v(310)*v(478)
       v(315)=-(v(312)*v(313))+v(311)*v(314)
       v(470)=d(2)*gp(4,i245)*v(315)
       v(321)=-(v(314)/v(315))
       v(335)=-(v(310)*v(321))
       v(327)=-(v(307)*v(321))
       v(322)=v(312)/v(315)
       v(338)=-(v(310)*v(322))
       v(329)=-(v(307)*v(322))
       v(323)=-(v(313)/v(315))
       v(336)=v(308)*v(323)
       v(331)=v(309)*v(323)
       v(324)=v(311)/v(315)
       v(339)=v(308)*v(324)
       v(333)=v(309)*v(324)
       v(325)=v(327)+v(336)
       v(326)=v(329)+v(339)
       v(328)=-v(327)+v(331)
       v(330)=-v(329)+v(333)
       v(332)=-v(331)+v(335)
       v(334)=-v(333)+v(338)
       v(337)=-v(335)-v(336)
       v(341)=v(281)*v(325)+v(282)*v(328)+v(283)*v(332)+v(284)*v(337)
       v(340)=-v(338)-v(339)
       v(342)=v(281)*v(326)+v(282)*v(330)+v(283)*v(334)+v(284)*v(340)
       p(1)=p(1)+(v(325)*v(341)+v(326)*v(342))*v(470)
       p(2)=p(2)+(v(328)*v(341)+v(330)*v(342))*v(470)
       p(3)=p(3)+(v(332)*v(341)+v(334)*v(342))*v(470)
       p(4)=p(4)+(v(337)*v(341)+v(340)*v(342))*v(470)
      ENDDO
      END

!******************* S U B R O U T I N E **********************
      SUBROUTINE elmt13_ISW09(v,d,xl,ul,c,p,ht,hp,gp,ngpo)
      IMPLICIT NONE
      include 'sms.h'
      INTEGER ngpo,i353
      DOUBLE PRECISION v(592),d(3),xl(2,4),ul(6,4),c(4,4),p
     &(4),ht(0),hp(0),gp(4,4)
      DO i353=1,ngpo
       c(1,1)=0d0+c(1,1)
       c(1,2)=0d0+c(1,2)
       c(1,3)=0d0+c(1,3)
       c(1,4)=0d0+c(1,4)
       c(2,1)=0d0+c(2,1)
       c(2,2)=0d0+c(2,2)
       c(2,3)=0d0+c(2,3)
       c(2,4)=0d0+c(2,4)
       c(3,1)=0d0+c(3,1)
       c(3,2)=0d0+c(3,2)
       c(3,3)=0d0+c(3,3)
       c(3,4)=0d0+c(3,4)
       c(4,1)=0d0+c(4,1)
       c(4,2)=0d0+c(4,2)
       c(4,3)=0d0+c(4,3)
       c(4,4)=0d0+c(4,4)
      ENDDO
      END
