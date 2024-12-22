      subroutine CoefLoadobs(BLH,hd,BB,nx,GF,dr,GRS,sgn,sk,kk)
!sgn=1高程异常mm,2扰动重力uGal,3地面重力uGal,4大地高mm,5正常高mm变化
!2015年8月1日,章传银
      implicit none
	integer::nx,sgn,sk(nx),kk
	real*8::BB(nx),hd(6),BLH(3),GF(8000,9),GRS(6),dr,gg,tmp,xx,yy,dd
 	real*8::dlat,dlon,ae,pi,gm,em,gr,RAD,vfn(9),dl,ds,rw,re,xs(10),rr,r1,tt
      real*8::sin2f,psi,rln(3),XYZ(3),BLH1(3),rln1(3),XYZ1(3),mdr
	integer::i,j,ni,nj,i0,j0,nlat,nlon,ki
!---------------------------------------------------------------------------
   	gg=6.67428d-11;gm=GRS(1); ae=GRS(2); rw=1.d3;em=gm/gg
   	gr=9.7803278d0;pi=datan(1.d0)*4.d0;RAD=pi/180.d0;BB=0.d0
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
	nlat=nint((hd(4)-hd(3))/hd(6));nlon=nint((hd(2)-hd(1))/hd(5))
      ni=nint(dr/ae/RAD/hd(6)+1.d0);nj=nint(dr/ae/RAD/hd(5)/dcos(rln(2)*RAD)+1.d0) !积分半径dr对应的地面格网数
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
      mdr=hd(5)*RAD*ae/4.d0;kk=0
	do i=i0-ni,i0+ni
	  if(i<1.or.i>nlat)goto 1001
        BLH1(1)=hd(3)+(dble(i)-0.5d0)*hd(6)
	  do j=j0-nj,j0+nj
	    if(j<1.or.j>nlon)goto 1002
	    BLH1(2)=hd(1)+(dble(j)-0.5d0)*hd(5)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          dlat=(rln1(2)-rln(2))*RAD;dlon=(rln1(3)-rln(3))*RAD
          xx=2.d0*dsin(dlat/2.d0)*r1
          yy=2.d0*dsin(dlon/2.d0)*r1*dcos((rln(2)+rln1(2))*RAD/2.d0)
          dl=dsqrt(xx**2+yy**2);if(dl>dr)goto 1002
          sin2f=dl/r1/2.d0
          ds=hd(5)*hd(6)*dcos(rln1(2)*RAD)*RAD**2*r1**2
          tt=1.d0-2.d0*sin2f**2; call IntpGrnF(GF,dl,vfn)
          dd=dl;if(dl<mdr)dd=mdr
          !直接影响               
          xs(1)=gg*rw*ds/dd/gr*1.d3
          xs(2)=gg*rw*ds/dd**3*(rr-r1*tt)*1.0e8
          !间接影响
          xs(1)=xs(1)+rw*ds*vfn(1)/dl*1.d3
          xs(2)=xs(2)+rw*ds*vfn(2)/dl*1.0e8
          xs(4)=rw*ds*vfn(7)/dl*1.d3
          xs(3)=xs(2)-xs(4)*0.3086d0
          xs(5)=xs(4)-xs(1)
          ki=(i-1)*nlon+j;kk=kk+1;sk(kk)=ki
          BB(ki)=xs(sgn)
1002      continue
        enddo
1001    continue
      enddo
 	return
      end
