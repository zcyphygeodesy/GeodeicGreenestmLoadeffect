!  GeodeicGreenestmLoadeffect.f90 
!
!  FUNCTIONS:
!  GeodeicGreenestmLoadeffect - Entry point of console application.
!
!****************************************************************************
      program GeodeicGreenestmLoadeffect
      implicit none
	character*800::dtmfl,loadgrfl,obstmsqdfl
	character*80000::line,headmn
	integer kln,sn,astat(15),tmlth
	integer knd,mode,itern,i,hepch,frow,kndrow,wghrow,nfar
	real*8::rec(8000),inp(12),lnth,dr,lvb,kwgh,GF(8000,9)
	integer::status=0
!---------------------------------------------------------------------
      !read load green functions for indirect load effect 负荷格林函数(间接影响)
      write(loadgrfl,*) "LoadGreen.txt"
      call LGrnFunc(loadgrfl,GF)
      !输入计算面高度格网文件名
      !Input the heterogeneous geodetic variation record time series file name.
      write(dtmfl,*)'dtm3m.dat'
      !输入大地测量站点记录时间序列文件名
      !Input the calculation surface height grid file name.
      write(obstmsqdfl,*)'heterobstm.txt'
      !输入站点平均间距lnth(m)和格林函数积分半径dr(m)
      !Input mean distance (m) between geodetic sites and Green's integral radius (m)
      lnth=20.d3;dr=150.d3
      !设置Laplace平滑算子权值与累计逼近次数
      !Set Laplace operator weight and cumulative approach times
      lvb=0.3d0;itern=2
      nfar=1!设置边缘效应抑制参数 Set edge effect suppression parameter
      !设置可调控监测量类型及其贡献率
      !Set type (1~5) of adjustable variation and contribution rate
      mode=3;kwgh=1.d0!mode=1~5;kwgh=1.d0 意味不调控
      !输入法方程解算方法 the method of the solution of normal equation
      !knd=1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square qr decomposition,
      !   =4 Minimum norm singular value decomposition
      knd=1!knd=1-LU分解,2-Cholesky分解,3-最小二乘QR分解,4-最小范数奇异值分解
      hepch=2!记录时序头文件首个采样时刻列序号 Column ordinal number of the first epoch time in header
      frow=7!记录时序中首次采样列序号 Column ordinal number of the first variation in record
      kndrow=6!监测量类型列序号 Column ordinal number of the variation type in record
      wghrow=5!监测量权值列序号 !Column ordinal number of the variation weight in record
      !将输入参数综合成向量inp(12) merge input parameters into inp(12)
      inp(1)=hepch;inp(2)=frow;inp(3)=kndrow;inp(4)=wghrow
      inp(5)=nfar;inp(6)=lnth;inp(7)=dr;inp(8)=lvb;inp(9)=kwgh
      inp(10)=knd;inp(11)=mode;inp(12)=itern
      open(unit=8,file=obstmsqdfl,status="old",iostat=status)
      if(status/=0)goto 902
      read(8,'(a)') headmn
      call PickReclong(headmn,kln,rec,sn)
      if(sn<hepch)then
         close(8);goto 902!缺采样历元时刻信息 Missing epoch time information
      endif
      close(8)
      tmlth=sn-hepch+1!监测量时序长度 length of geodetic variation record time series
      write(*, *)"    Begin compulation......"
      do i=frow,frow+tmlth-1 !i=当前监测量列序号 column ordinal number of the current variations in record
        call Loadeffectestmgreen(dtmfl,obstmsqdfl,GF,i,inp)
        write(*, 104)nint(rec(hepch+i-frow))
      enddo
902   continue
      pause
104   format('     Load EWH and load effects computed at Epoch time',I12)
      end
