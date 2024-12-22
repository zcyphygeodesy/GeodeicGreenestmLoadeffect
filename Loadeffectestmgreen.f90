      subroutine Loadeffectestmgreen(dtmfl,obstmsqdfl,GF,dtrow,inp)
      USE f95_precision, ONLY: WP => DP
      USE lapack95
      implicit none
	character*800::dtmfl,obstmsqdfl,viewfl,outfl,rntfl
	character*80000::line,headmn,nmstr(5),stm,pstr
      character(len=25)::str(8000)
	integer dtrow,i,j,k,n,m,nn,kk,sn,kln,astat(15),nx,ni,mi,nfar,itern,mxk,sgn,kp,knd,kobs,mode
 	real*8 GRS(6),inp(12),hd(6),hd0(6),wgh,lon,lat,dr,RAD,rec(8000),st(5,4,5),wght(5)
 	real*8 GF(8000,9),BLH(3),tdn(24),tsr(4),ldn(14),lvb,bf(8),tmp,lnth,ab(5),sta(5,4),Gauss2D
	integer nlon,nlat,nlon0,nlat0,nk,maxn,hepch,frow,kndrow,wghrow,tmlth,obsn(5)!
	real*8,allocatable::BPB(:,:,:),BPB0(:,:),BPL(:,:),s(:),xx(:),ewh(:,:),ewh1(:,:)
	real*8,allocatable::obs(:,:,:),APA(:,:),APL(:),tm(:),chs(:,:,:),rst(:,:,:)
	integer,allocatable::sk(:)
	character(len=25),allocatable::strnm(:,:)
	integer::status=0
!kobs=1高程异常mm,2扰动重力uGal,3地面重力uGal,4大地高mm,5正常高mm变化
!输出ewh,geoid,terrgrav,gravdist,grndtilt,vertdefl,horzdisp,elliphgt,orthohgt,gradient,horzgrad
!---------------------------------------------------------------------
      write(nmstr(1),'(a20,I3)')'Height anomaly:',1  !mm
      write(nmstr(2),'(a20,I3)')'Gravity disturbance:',2  !uGal
      write(nmstr(3),'(a20,I3)')'Ground gravity:',3   !uGal
      write(nmstr(4),'(a20,I3)')'Ellipsoidal height:',4   !mm
      write(nmstr(5),'(a20,I3)')'Normal Height:',5  !mm
    	RAD=datan(1.d0)/45.d0;BLH(3)=0.d0;wght(1:5)=1.d0
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0
      hepch=nint(inp(1));frow=nint(inp(2));kndrow=nint(inp(3));wghrow=nint(inp(4))
      nfar=nint(inp(5));lnth=inp(6);dr=inp(7);lvb=inp(8);knd=nint(inp(10))
      mode=nint(inp(11));itern=nint(inp(12));wght(mode)=inp(9)
      open(unit=8,file=dtmfl,status="old",iostat=status)
      if(status/=0)goto 902
	read(8,*,end=905)(hd0(i),i=1,6)
905   close(8)
	nlat0=nint((hd0(4)-hd0(3))/hd0(6));nlon0=nint((hd0(2)-hd0(1))/hd0(5))
	hd0(5)=(hd0(2)-hd0(1))/dble(nlon0);hd0(6)=(hd0(4)-hd0(3))/dble(nlat0)
      hd(1:4)=hd0(1:4);hd(5)=lnth/GRS(2)/RAD/2.d0
      nlat=nint((hd(4)-hd(3))/hd(5));nlon=nint((hd(2)-hd(1))/hd(5))
      hd(5)=(hd(4)-hd(3))/dble(nlat);hd(6)=(hd(2)-hd(1))/dble(nlon)
      nx=nlat*nlon; nn=0 !nn为有效地面测站点数
      open(unit=8,file=obstmsqdfl,status="old",iostat=status)
      if(status/=0)goto 902
      read(8,'(a)') headmn
      call PickReclong(headmn,kln,rec,sn)
      nn=0!有效观测数
      if(sn<hepch)then
         close(8);goto 902!缺采样历元时刻信息
      endif
      tmlth=sn-hepch+1!采样时序长度
      do while(.not.eof(8))
        read(8,'(a)') line
        call PickReclong(line,kln,rec,sn)
        if(rec(2)<hd(1).or.rec(2)>hd(2).or.rec(3)<hd(3).or.rec(3)>hd(4))goto 505
        if(sn<kndrow.or.sn<wghrow.or.sn<frow.or.sn<dtrow)goto 505
        if(sn-frow+1<tmlth)goto 505!缺监测量采样时序
        nn=nn+1
 505    continue
	enddo
      close(8)
      if(nn<6)goto 902 !监测量太少!中止计算！
 	allocate(BPB(nx,nx,5), stat=astat(1))
 	allocate(BPL(nx,5), stat=astat(2))
 	allocate(s(nx), stat=astat(3))
 	allocate(tm(tmlth), stat=astat(4))
 	allocate(sk(nx), stat=astat(5))
 	allocate(ewh(nlat,nlon), stat=astat(6))
 	allocate(ewh1(nlat,nlon), stat=astat(7))
 	allocate(BPB0(nx,nx), stat=astat(8))
 	allocate(APA(nx,nx), stat=astat(9))
 	allocate(APL(nx), stat=astat(10))
 	allocate(obs(nn,5,5), stat=astat(11))!5种观测场元-纬度经度，观测量，权值，原观测量
 	allocate(chs(nn,5,5), stat=astat(12))
 	allocate(strnm(nn,5), stat=astat(13))
 	allocate(rst(nlat0,nlon0,15), stat=astat(14))
	if (sum(astat(1:14)) /= 0)goto 902 !数字高程模型分辨率太大！计算中断！
      open(unit=8,file=obstmsqdfl,status="old",iostat=status)
      read(8,'(a)') headmn
      call PickReclong(headmn,kln,rec,sn)
      tm(1:tmlth)=rec(hepch:sn)
      nn=0;obsn=0
      do while(.not.eof(8))
        read(8,'(a)') line
        call PickReclong(line,kln,rec,sn)
        if(rec(2)<hd(1).or.rec(2)>hd(2).or.rec(3)<hd(3).or.rec(3)>hd(4))goto 506
        if(sn<kndrow.or.sn<wghrow.or.sn<frow.or.sn<dtrow)goto 506
        if(sn-frow+1<tmlth)goto 506!缺监测量采样时序
        if(rec(kndrow)<0.5d0.or.rec(kndrow)>5.5d0)goto 506!监测量类型超限
        nn=nn+1;nk=nint(rec(kndrow));obsn(nk)=obsn(nk)+1
        obs(obsn(nk),1,nk)=rec(3);obs(obsn(nk),2,nk)=rec(2)
        obs(obsn(nk),3,nk)=rec(dtrow);obs(obsn(nk),4,nk)=rec(wghrow)
        obs(obsn(nk),5,nk)=rec(dtrow)!5原观测量
        call PickRstrlg(line,kln,str,sn)
        strnm(obsn(nk),nk)=str(1)
 506    continue
	enddo
      close(8)
      kp=0;ewh1=0.d0;BLH(3)=0.d0
1011  BPB=0.d0;BPL=0.d0
      do kobs=1,5
        do k=1,obsn(kobs)
          tmp=obs(k,3,kobs)
	    if(tmp>9.d3)goto 2001
          BLH(1:2)=obs(k,1:2,kobs);wgh=obs(k,4,kobs);sk=0
          call CoefLoadobs(BLH,hd,s,nx,GF,dr,GRS,kobs,sk,mxk)!sk(mxk)存储是s(sk(mxk))
          do n=1,mxk
             ni=sk(n);BPL(ni,kobs)=BPL(ni,kobs)+s(ni)*tmp*wgh
             if(kp>0)goto 2002
             do m=1,mxk
               mi=sk(m);BPB(ni,mi,kobs)=BPB(ni,mi,kobs)+s(ni)*s(mi)*wgh
             enddo
2002         continue
          enddo
2001      continue
        enddo
        if(obsn(kobs)>0.and.kp==0)call Stat1d(obs(1:obsn(kobs),3,kobs),obsn(kobs),sta(kobs,1:4))
      enddo
      if(kp>0)then
        APA=BPB0; goto 2003
      endif
      APA=0.d0;ab=0.d0
      do k=1,5!计算法方程对角线标准差
        if(obsn(k)<1)goto 3131
        s=0.d0
        do i=1,nx
           s(i)=BPB(i,i,k)
        enddo
        call Stat1d(s(1:nx),nx,tsr)
        ab(k)=3.d0*tsr(2)+(maxval(s(1:nx))-minval(s(1:nx)))*0.1d0
3131    continue
      enddo
      do k=1,5
        if(obsn(k)<1)goto 3232
        tmp=wght(k)/ab(k)!法方程对角线标准差法组合
        do i=1,nx
          do j=1,nx
            APA(i,j)=APA(i,j)+tmp*BPB(i,j,k)
	    enddo
        enddo
3232    continue
      enddo
      !增加Laplace平滑算子，配权lvb
1601  call Laplacesmth(APA,s,nlat,nlon,lvb)
      !增加四周远区零效应
      tmp=0.d0
      do i=1,nx
        tmp=tmp+APA(i,i)**2/dble(nx)
      enddo
      tmp=dsqrt(tmp)
      call Farzonexx(APA,nlat,nlon,nfar,tmp)
2003  APL=0.d0
      do k=1,5
        if(obsn(k)<1)goto 3233
        tmp=wght(k)/ab(k)!法方程对角线标准差法组合
        do i=1,nx
          APL(i)=APL(i)+tmp*BPL(i,k)
        enddo
3233    continue
      enddo
      BPB0=APA;s=0.d0
      !knd=1-LU分解,2-Cholesky分解,3-最小二乘QR分解,4-最小范数奇异值分解
      call Equsolve(APA,s,nx,APL,knd,bf)
	do i=1,nlat
         do j=1,nlon
            k=(i-1)*nlon+j;ewh(i,j)=s(k)
         enddo
      enddo
      do kobs=1,5
        do k=1,obsn(kobs)
          tmp=obs(k,3,kobs)
	    if(tmp>9.d3)goto 2004
          BLH(1:2)=obs(k,1:2,kobs)
          call rntGreenintgr(BLH,ewh,hd,nlat,nlon,GF,ldn,GRS,dr)
          tdn(1:3)=ldn(1:3);tdn(4)=ldn(10);tdn(5)=tdn(4)-tdn(1)
          chs(k,kobs,kp+1)=obs(k,3,kobs)-tdn(kobs);obs(k,3,kobs)=chs(k,kobs,kp+1)
2004      continue
        enddo
        if(obsn(kobs)>0)call Stat1d(chs(1:obsn(kobs),kobs,kp+1),obsn(kobs),st(kobs,1:4,kp+1))
      enddo
      kp=kp+1;ewh1=ewh1+ewh
      write(stm,*)nint(tm(dtrow-frow+1))
 	write(rntfl,*) "testdata\rnt",trim(AdjustL(stm)),".txt"
      open(unit=10,file=rntfl,status="replace")
      do kobs=1,5
        if(obsn(kobs)>1)then
           write(10,'(a24)')trim(nmstr(kobs))
           write(10,'(I8,4F10.4)')0,(sta(kobs,i),i=1,4)
           do i=1,kp
              write(10,'(I8,4F10.4)')i, st(kobs,1:4,i)
           enddo
        endif
      enddo
      do kobs=1,5
        do k=1,obsn(kobs)
          BLH(1:2)=obs(k,1:2,kobs);write(str(1),'(a)')strnm(k,kobs)
	    write(10,'(a12,2F10.4,F8.2,I3,80F12.4)')trim(str(1)),BLH(2),BLH(1),obs(k,4,kobs),kobs,obs(k,5,kobs),chs(k,kobs,1:kp)
        enddo
      enddo
      close(10)
      if(kp<itern+1) goto 1011
	do i=1,nlat0
        BLH(1)=hd0(3)+(real(i)-0.5d0)*hd0(6)
        do j=1,nlon0
           BLH(2)=hd0(1)+(real(j)-0.5d0)*hd0(5)
           rst(i,j,1)=Gauss2D(BLH(2),BLH(1),ewh1,nlat,nlon,hd)
        enddo
	enddo
	do i=1,nlat0
        BLH(1)=hd0(3)+(real(i)-0.5d0)*hd0(6)
        do j=1,nlon0
           BLH(2)=hd0(1)+(real(j)-0.5d0)*hd0(5)
           call rntGreenintgr(BLH,rst(:,:,1),hd0,nlat0,nlon0,GF,ldn,GRS,dr)!ldn(14)
           rst(i,j,2:15)=ldn(1:14)
        enddo
	enddo
 	write(outfl,*) "testdata\ewh",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,1)*1.d2,j=1,nlon0)!等效水高m->cm
 	enddo
      close(10)
 	write(outfl,*) "testdata\Greengeoid",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,2),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*)"testdata\Greenterrgrav",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,3),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*) "testdata\Greengravdist",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,4),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*) "testdata\Greengrndtilt",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,5),j=1,nlon0)
 	enddo
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,6),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*) "testdata\Greenvertdefl",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,7),j=1,nlon0)
 	enddo
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,8),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*)"testdata\Greenhorzdisp",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,9),j=1,nlon0)
 	enddo
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,10),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*)"testdata\Greenelliphgt",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,11),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*) "testdata\Greenorthohgt",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15F12.4)')(rst(i,j,12),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*)"testdata\Greengradient",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15ES14.5)')(rst(i,j,13),j=1,nlon0)
 	enddo
      close(10)
 	write(outfl,*)"testdata\Greenhorzgrad",trim(AdjustL(stm)),".dat"
      open(unit=10,file=outfl,status="replace")
      write(10,101)(hd0(i),i=1,6),nint(tm(dtrow-frow+1))
	do i=1,nlat0
	  write(10,'(15ES14.5)')(rst(i,j,14),j=1,nlon0)
 	enddo
	do i=1,nlat0
	  write(10,'(15ES14.5)')(rst(i,j,15),j=1,nlon0)
 	enddo
      close(10)
901 	deallocate(ewh,ewh1,BPB,BPB0,BPL,s,tm,obs,sk,chs,APA,APL,strnm,rst)
902   continue
101   format(4F12.6,2ES16.8,I14)
      end
!
!***********************************************************************
!
      subroutine Laplacesmth(BPB,s,nlat,nlon,wgh)
      !Laplace平滑算子,BPB法方程-修订，wgh权,s(nlat*nlon)提供空间
      !spatial Laplace filter constraints supplemented on unknown parameters
!-------------------------------------------------------------
      implicit none
	integer::i,j,k,nlat,nlon,n,m,lp(5)
	real*8::s(nlat*nlon),BPB(nlat*nlon,nlat*nlon),wgh
!-----------------------------------------------------------------
	do i=1,nlat
         do j=1,nlon
           s=0.d0;k=0!Laplace算子
           if(i-1>0)then
              s((i-2)*nlon+j)=1.d0;k=k+1;lp(k)=(i-2)*nlon+j
           endif
           if(i<nlat)then
              s(i*nlon+j)=1.d0;k=k+1;lp(k)=i*nlon+j
           endif
           if(j-1>0)then
              s((i-1)*nlon+j-1)=1.d0;k=k+1;lp(k)=(i-1)*nlon+j-1
           endif
           if(j<nlon)then
              s((i-1)*nlon+j+1)=1.d0;k=k+1;lp(k)=(i-1)*nlon+j+1
           endif
           s((i-1)*nlon+j)=-dble(k);k=k+1;lp(k)=(i-1)*nlon+j
           do n=1,k
             do m=1,k
               BPB(lp(n),lp(m))=BPB(lp(n),lp(m))+s(lp(n))*s(lp(m))*wgh
             enddo
           enddo
         enddo
      enddo
      end
!
!************************************************************************
!
      subroutine Farzonexx(BPB,nlat,nlon,n,wgh)
      !四周远区零效应，BPB法方程-修订，wgh权
      !far-zone zero constraint supplemented on unknown parameters at grid edge
!-------------------------------------------------------------
      implicit none
	integer::i,j,k,nlat,nlon,n
	real*8::BPB(nlat*nlon,nlat*nlon),wgh
!-----------------------------------------------------------------
	do i=1,n
        do j=1,nlon
          k=(i-1)*nlon+j;BPB(k,k)=BPB(k,k)+wgh
        enddo
      enddo
	do i=nlat-n+1,nlat
        do j=1,nlon
          k=(i-1)*nlon+j;BPB(k,k)=BPB(k,k)+wgh
        enddo
      enddo
	do j=1,n
        do i=1,nlat
          k=(i-1)*nlon+j;BPB(k,k)=BPB(k,k)+wgh
        enddo
      enddo
	do j=nlon-n+1,nlon
        do i=1,nlat
          k=(i-1)*nlon+j;BPB(k,k)=BPB(k,k)+wgh
        enddo
      enddo
      end
