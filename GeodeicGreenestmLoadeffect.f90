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
      !read load green functions for indirect load effect ���ɸ��ֺ���(���Ӱ��)
      write(loadgrfl,*) "LoadGreen.txt"
      call LGrnFunc(loadgrfl,GF)
      !���������߶ȸ����ļ���
      !Input the heterogeneous geodetic variation record time series file name.
      write(dtmfl,*)'dtm3m.dat'
      !�����ز���վ���¼ʱ�������ļ���
      !Input the calculation surface height grid file name.
      write(obstmsqdfl,*)'heterobstm.txt'
      !����վ��ƽ�����lnth(m)�͸��ֺ������ְ뾶dr(m)
      !Input mean distance (m) between geodetic sites and Green's integral radius (m)
      lnth=20.d3;dr=150.d3
      !����Laplaceƽ������Ȩֵ���ۼƱƽ�����
      !Set Laplace operator weight and cumulative approach times
      lvb=0.3d0;itern=2
      nfar=1!���ñ�ԵЧӦ���Ʋ��� Set edge effect suppression parameter
      !���ÿɵ��ؼ�������ͼ��乱����
      !Set type (1~5) of adjustable variation and contribution rate
      mode=3;kwgh=1.d0!mode=1~5;kwgh=1.d0 ��ζ������
      !���뷨���̽��㷽�� the method of the solution of normal equation
      !knd=1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square qr decomposition,
      !   =4 Minimum norm singular value decomposition
      knd=1!knd=1-LU�ֽ�,2-Cholesky�ֽ�,3-��С����QR�ֽ�,4-��С��������ֵ�ֽ�
      hepch=2!��¼ʱ��ͷ�ļ��׸�����ʱ������� Column ordinal number of the first epoch time in header
      frow=7!��¼ʱ�����״β�������� Column ordinal number of the first variation in record
      kndrow=6!�������������� Column ordinal number of the variation type in record
      wghrow=5!�����Ȩֵ����� !Column ordinal number of the variation weight in record
      !����������ۺϳ�����inp(12) merge input parameters into inp(12)
      inp(1)=hepch;inp(2)=frow;inp(3)=kndrow;inp(4)=wghrow
      inp(5)=nfar;inp(6)=lnth;inp(7)=dr;inp(8)=lvb;inp(9)=kwgh
      inp(10)=knd;inp(11)=mode;inp(12)=itern
      open(unit=8,file=obstmsqdfl,status="old",iostat=status)
      if(status/=0)goto 902
      read(8,'(a)') headmn
      call PickReclong(headmn,kln,rec,sn)
      if(sn<hepch)then
         close(8);goto 902!ȱ������Ԫʱ����Ϣ Missing epoch time information
      endif
      close(8)
      tmlth=sn-hepch+1!�����ʱ�򳤶� length of geodetic variation record time series
      write(*, *)"    Begin compulation......"
      do i=frow,frow+tmlth-1 !i=��ǰ���������� column ordinal number of the current variations in record
        call Loadeffectestmgreen(dtmfl,obstmsqdfl,GF,i,inp)
        write(*, 104)nint(rec(hepch+i-frow))
      enddo
902   continue
      pause
104   format('     Load EWH and load effects computed at Epoch time',I12)
      end
