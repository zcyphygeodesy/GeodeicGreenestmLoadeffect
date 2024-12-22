      subroutine rntGreenintgr(BLH,ewh,hd,nlat,nlon,GF,ldn,GRS,dr)
      !dr积分半径m
      implicit none
	integer::nlat,nlon,i,j,ni,nj,i0,j0
	real*8::ewh(nlat,nlon),BLH(3),hd(6),direct(10),tdn(8),ldn(14),GRS(6),GF(8000,9),indrct(14)
 	real*8::pi,gg,gr,RAD,L1,ds,rw,L0,mdr,dv,dwh,dr,rst(4),hp,kp,tmp,dl,dv0
      real*8::rln(3),XYZ(3),BLH1(3),vfn(9),rln1(3),XYZ1(3),sin2f,tt,NFD(5),dfi0
      real*8::BLH0(3),rln0(3),XYZ0(3),BLH2(3),rln2(3),XYZ2(3),rr,r0,r1
	real*8 rlon,rlat,rlon1,rlat1,dlon,dlat,cosa,sina,cos2f,sinf,cosf,ctgf
!---------------------------------------------------------------------------
      hp=1.d0-BLH(3)/4.43d4;if(hp<0.d0)hp=1.d-12;kp=dexp(dlog(hp)*5.256d0)
      kp=1.d0-2.d0*kp;if(dabs(GRS(6)+1.d0)>1.d-8)kp=1.d0!GRS(6)=-1.d0大气压负荷影响
      gg=6.67428d-11;rw=1.d3;pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      BLH0=BLH;BLH0(3)=0.d0;ldn=0.d0;tdn(1:14)=0.d0;indrct=0.d0;direct=0.d0
      call BLH_RLAT(GRS,BLH,rln);call BLH_XYZ(GRS,BLH,XYZ)
      rr=rln(1);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      rlat=rln(2)*RAD;rlon=rln(3)*RAD
      call BLH_XYZ(GRS,BLH0,XYZ0)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      mdr=r0*hd(5)*RAD*dcos(rlat)/2.d0 !奇异点判断
      ni=nint(dr/r0/RAD/hd(6)+1.d0)
      nj=nint(dr/r0/RAD/hd(5)/dcos(rlat)+1.d0) !积分半径dr对应的地面格网数
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
	do i=i0-ni,i0+ni
	  if(i<1.or.i>nlat)goto 9100
        BLH1(1)=hd(3)+(dble(i)-0.5d0)*hd(6)
	  do j=j0-nj,j0+nj
	    if(j<1.or.j>nlon)goto 9101
          if(dabs(ewh(i,j))<1.d-8)goto 9101
	    BLH1(2)=hd(1)+(dble(j)-0.5d0)*hd(5)
          BLH1(3)=0.d0;call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          if(L0>dr)goto 9101
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          sin2f=L0/r0/2.d0;tt=1.d0-2.d0*sin2f**2
          cos2f=dsqrt(1.d0-sin2f**2);sinf=2.d0*sin2f*cos2f
          cosf=dsqrt(1.d0-sinf**2)
          rlat1=rln1(2)*RAD;rlon1=rln1(3)*RAD
          ds=hd(5)*hd(6)*RAD**2*dcos(rlat1)*r1**2
          dlat=rlat1-rlat;dlon=rlon1-rlon
          dwh=gg*rw*ewh(i,j);dv=dwh*ds
          dl=L0;call IntpGrnF(GF,dl,vfn)
          if(L1<mdr)then !奇异点处理
            call sgnGrnintigral(BLH,ewh,hd,nlat,nlon,GF,tdn,indrct,GRS,dr)
          else
            dl=L0;call IntpGrnF(GF,dl,vfn)
	      cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(dlon))/sinf
	      sina=dcos(rlat1)*dsin(dlon)/sinf;ctgf=cosf/sinf
            tdn(1)=tdn(1)+dv/L1
            tdn(2)=tdn(2)+dv/L1**3*(rr-r1*tt)
            tdn(3)=tdn(3)-dv*r1/L1**3*cosa*sinf
            tdn(4)=tdn(4)-dv*r1/L1**3*sina*sinf*dcos(rlat)
            tdn(5)=tdn(5)-dv*(1.d0/L1**3-3.d0*(rr-r1*cosf)**2/L1**5)
            tmp=dv*r1/rr*(cosf/L1**3-3.d0*rr*r1*sinf**2/L1**5)
            tdn(6)=tdn(6)-tmp*ctgf*(1.d0-cosa**2)
            tdn(7)=tdn(7)+tmp*((1.d0-dcos(rlat)**2*sina**2)*ctgf-dcos(rlat1)*dcos(rlat)/sinf)
            dv0=rw*ewh(i,j)*ds
            indrct(10)=indrct(10)+dv0*vfn(7)/dl
            indrct(1)=indrct(1)+dv0*vfn(1)/dl
            indrct(3)=indrct(3)+dv0*vfn(3)/dl
            indrct(12)=indrct(12)+dv0*vfn(8)/dl
            indrct(4)=indrct(4)+dv0/dl*vfn(4)*cosa;indrct(5)=indrct(5)+dv0/dl*vfn(4)*sina*dcos(rlat1)
            indrct(6)=indrct(6)+dv0/dl*vfn(5)*cosa;indrct(7)=indrct(7)+dv0/dl*vfn(5)*sina*dcos(rlat1)
            indrct(8)=indrct(8)-dv0/dl*vfn(6)*sina*dcos(rlat);indrct(9)=indrct(9)-dv0/dl*vfn(6)*cosa
            indrct(13)=indrct(13)+dv0/dl*vfn(9)*ctgf*(1.d0-cosa**2)
            indrct(14)=indrct(14)-dv0/dl*vfn(9)*((1.d0-dcos(rlat)**2*sina**2)*ctgf-dcos(rlat1)*dcos(rlat)/sinf)
          endif
9101      continue
        enddo
9100    continue
      enddo
	direct(1)=tdn(1)/gr*1.d3
	direct(2)=tdn(2)*1.0e8*kp
      direct(3)=direct(2)
	direct(4)=tdn(3)/gr*36.d5/RAD
	direct(5)=tdn(4)/gr/dcos(rlat)*36.d5/RAD
      direct(6:7)=direct(4:5)
      direct(8)=-tdn(5)*1.d12*kp!重力梯度直接影响径向mE
      direct(9:10)=tdn(6:7)*1.d9!水平梯度直接影响E
      direct(10)=-direct(10)!西向
	indrct(1)=indrct(1)*1.d3
	indrct(3)=indrct(3)*1.0e8
	indrct(4)=indrct(4)/RAD*36.d5
	indrct(5)=indrct(5)/RAD*36.d5
	indrct(6)=indrct(6)/RAD*36.d5
	indrct(7)=indrct(7)/RAD*36.d5
	indrct(8:10)=indrct(8:10)*1.d3
	indrct(12)=-indrct(12)*1.d12!重力梯度间接影响径向mE
	indrct(13:14)=indrct(13:14)*1.d9!水平梯度间接影响E
      indrct(14)=-indrct(14)!西向
	indrct(2)=indrct(3)-indrct(10)*0.3086d0
      ldn(1:7)=direct(1:7)+indrct(1:7);ldn(8:10)=indrct(8:10)
      ldn(11)=ldn(10)-ldn(1);ldn(12:14)=direct(8:10)+indrct(12:14)
 	return
      end
!
!************************************************************************************
! 
      subroutine sgnGrnintigral(BLH,ewh,hd,nlat,nlon,GF,tdn,indrct,GRS,dr)
      !dr积分半径m
      implicit none
	integer::nlat,nlon,i,j,ni,nj,i0,j0
	real*8::ewh(nlat,nlon),BLH(3),hd(6),tdn(8),GRS(6),GF(8000,9),indrct(14)
 	real*8::pi,RAD,L1,ds,rw,L0,mdr,dv,dwh,dr,rst(4),hp,kp,tmp,dl,dv0,BLHz(3)
      real*8::rln(3),XYZ(3),BLH1(3),vfn(9),rln1(3),XYZ1(3),sin2f,tt,dfi0,gg,ewh0
      real*8::BLH0(3),rln0(3),XYZ0(3),BLH2(3),rln2(3),XYZ2(3),rr,r0,r1
	real*8 rlon,rlat,rlon1,rlat1,dlon,dlat,cosa,sina,cos2f,sinf,cosf,ctgf
!---------------------------------------------------------------------------
      hp=1.d0-BLH(3)/4.43d4;if(hp<0.d0)hp=1.d-12;kp=dexp(dlog(hp)*5.256d0)
      kp=1.d0-2.d0*kp;if(dabs(GRS(6)+1.d0)>1.d-8)kp=1.d0!GRS(6)=-1.d0大气压负荷影响
      gg=6.67428d-11;rw=1.d3;pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      BLH0=BLH;BLH0(3)=0.d0
      call BLH_RLAT(GRS,BLH,rln);call BLH_XYZ(GRS,BLH,XYZ)
      rr=rln(1)
      rlat=rln(2)*RAD;rlon=rln(3)*RAD
      call BLH_XYZ(GRS,BLH0,XYZ0)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      mdr=r0*hd(5)*RAD*dcos(rlat)/8.d0 !奇异点判断
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
      !计算格网左下角经纬度
      ewh0=ewh(i0,j0);BLHz(1)=hd(3)+dble(i0)*hd(6);BLHz(2)=hd(1)+dble(j0)*hd(5)
	do i=1,4
        BLH1(1)=BLHz(1)+(dble(i)-0.5d0)*hd(6)/4.d0
	  do j=1,4
	    BLH1(2)=BLHz(2)++(dble(j)-0.5d0)*hd(5)/4.d0
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          sin2f=L0/r0/2.d0;tt=1.d0-2.d0*sin2f**2
          cos2f=dsqrt(1.d0-sin2f**2);sinf=2.d0*sin2f*cos2f
          cosf=dsqrt(1.d0-sinf**2)
          rlat1=rln1(2)*RAD;rlon1=rln1(3)*RAD
          ds=hd(5)*hd(6)*RAD**2*dcos(rlat1)*r1**2/16.d0
          dlat=rlat1-rlat;dlon=rlon1-rlon
          dwh=gg*rw*ewh0;dv=dwh*ds
          if(L1<mdr)then
            L1=mdr;L0=mdr/2.d0;sinf=mdr/r0/2.d0
          endif
          dl=L0;call IntpGrnF(GF,dl,vfn)
	    cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(dlon))/sinf
	    sina=dcos(rlat1)*dsin(dlon)/sinf;ctgf=cosf/sinf
          tdn(1)=tdn(1)+dv/L1
          tdn(2)=tdn(2)+dv/L1**3*(rr-r1*tt)
          tdn(3)=tdn(3)-dv*r1/L1**3*cosa*sinf
          tdn(4)=tdn(4)-dv*r1/L1**3*sina*sinf*dcos(rlat)
          tdn(5)=tdn(5)-dv*(1.d0/L1**3-3.d0*(rr-r1*cosf)**2/L1**5)
          tmp=dv*r1/rr*(cosf/L1**3-3.d0*rr*r1*sinf**2/L1**5)
          tdn(6)=tdn(6)-tmp*ctgf*(1.d0-cosa**2)
          tdn(7)=tdn(7)+tmp*((1.d0-dcos(rlat)**2*sina**2)*ctgf-dcos(rlat1)*dcos(rlat)/sinf)
          dv0=rw*ewh(i,j)*ds
          indrct(10)=indrct(10)+dv0*vfn(7)/dl
          indrct(1)=indrct(1)+dv0*vfn(1)/dl
          indrct(3)=indrct(3)+dv0*vfn(3)/dl
          indrct(12)=indrct(12)+dv0*vfn(8)/dl
          indrct(4)=indrct(4)+dv0/dl*vfn(4)*cosa;indrct(5)=indrct(5)+dv0/dl*vfn(4)*sina*dcos(rlat1)
          indrct(6)=indrct(6)+dv0/dl*vfn(5)*cosa;indrct(7)=indrct(7)+dv0/dl*vfn(5)*sina*dcos(rlat1)
          indrct(8)=indrct(8)-dv0/dl*vfn(6)*sina*dcos(rlat);indrct(9)=indrct(9)-dv0/dl*vfn(6)*cosa
          indrct(13)=indrct(13)+dv0/dl*vfn(9)*ctgf*(1.d0-cosa**2)
          indrct(14)=indrct(14)-dv0/dl*vfn(9)*((1.d0-dcos(rlat)**2*sina**2)*ctgf-dcos(rlat1)*dcos(rlat)/sinf)
        enddo
      enddo
 	return
      end
