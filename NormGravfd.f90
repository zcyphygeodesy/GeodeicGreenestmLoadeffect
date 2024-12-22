      subroutine GNormalfd(BLH,NFD,GRS)
!MS$ATTRIBUTES DLLEXPORT::GNORMALFD
      !NFD(5)正常重力位,正常重力,正常重力梯度,正常重力线方向,正常梯度方向
      implicit none
 	real*8::BLH(3),NFD(5),val(5),rln(3),pn(40),dp1(40),dp2(40)
      integer::i,j,n
	real*8::GRS(6),djn(80),pi,RAD,gm,ae,ww,rr,t,u,vdn(5),atr
!---------------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      call normdjn(GRS,djn)
	gm=GRS(1);ae=GRS(2);ww=GRS(4);vdn=0.d0
      call BLH_RLAT(GRS,BLH,rln)
      rr=rln(1);t=dsin(rln(2)*RAD);u=dcos(rln(2)*RAD)
      call LegPn_dt2(pn,dp1,dp2,40,t)
      do n=1,20
        atr=dexp(dble(2*n)*dlog(ae/rr))*djn(2*n)
        vdn(1)=vdn(1)+atr*pn(2*n)
        vdn(2)=vdn(2)+atr*pn(2*n)*dble(2*n+1)
        vdn(3)=vdn(3)+atr*dp1(2*n)
        vdn(4)=vdn(4)+atr*pn(2*n)*dble((2*n+1)*(n+1))
        vdn(5)=vdn(5)+atr*dp2(2*n)
      enddo
      val(1)=gm/rr*(1.d0-vdn(1))+(ww*rr*u)**2/2.d0
      val(2)=-gm/rr**2*(1.d0-vdn(2))+ww**2*rr*u**2 !径向正常重力
      val(3)=-gm/rr**2*vdn(3)+ww**2*rr*u*t         !南向正常重力
      val(4)=-2.d0*gm/rr**3*(1.d0-vdn(4))+(ww*u)**2
      val(5)=-gm/rr**3*vdn(5)+ww**2*(2.d0*t**2-1.d0)
      NFD(4)=datan(val(3)/val(2))/RAD*36.d2
      NFD(5)=datan(val(5)/val(4))/RAD*36.d2
      val(2)=dsqrt(val(2)**2+val(3)**2)
      val(4)=dsqrt(val(4)**2+val(5)**2)
      NFD(1)=val(1);NFD(2)=val(2);NFD(3)=val(4)
      end
!
!============================================================================
!
      subroutine normdjn(GRS,djn)
!MS$ATTRIBUTES DLLEXPORT::NORMDJN
!输入GRS－gm,ae,j2,omega,1/f
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     djn  array of even zonals of normal ellipsoid
!c     ex. djn(2)=j2=GRS(3)
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)
      real*8 GRS(6),djn(80)
	ae=GRS(2)
	ff=GRS(5)
	do 101 n2=1,80,1
	djn(n2)=0.d0
  101	continue
	djn(2)=GRS(3)
      esq=(2.d0-ff)*ff!e**2
!c     compute normal even zonals from 4 to 20
      do 200 n2=4,80,2
        n=n2/2
        djn(n2)=(-1.d0)**(n+1)*3.d0*esq**n/(n2+1.d0)/(n2+3.d0)
     >        *(1.d0-dble(n)+5.d0*dble(n)*djn(2)/esq)
  200 continue
      return
      end
!
!******************************************************************************
!
      subroutine LegPn_dt2(pn,dp1,dp2,n,t)
      !计算Pn(t)及其对θ一、二阶导数t=cosθ
	implicit none
	integer::n,k
	real*8::pn(n),dp1(n),dp2(n),t,u
!---------------------------------------------------------------------------
      u=dsqrt(1-t**2)
      pn(1)=t;pn(2)=1.5d0*t**2-0.5d0
      dp1(1)=-u;dp1(2)=-3.d0*t*u
      dp2(1)=-t;dp2(2)=3.d0*(1.d0-2.d0*t**2)
      do k=3,n
        pn(k)=dble(2*k-1)/dble(k)*t*pn(k-1)-dble(k-1)*pn(k-2)/dble(k)
        dp1(k)=dble(2*k-1)*(t*dp1(k-1)-u*pn(k-1))-dble(k-1)*dp1(k-2)
        dp1(k)=dp1(k)/dble(k)
        dp2(k)=(t*dp2(k-1)-2.d0*u*dp1(k-1)-t*pn(k-1))*dble(2*k-1)
        dp2(k)=(dp2(k)-dble(k-1)*dp2(k-2))/dble(k)
      enddo
      end
