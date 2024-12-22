       subroutine LGrnFunc(loadgrfl,GF)
      !读取负荷格林函数(1:1000),dl=100m
      !GF(:,7)-高程异常,地面重力,扰动重力,地倾斜psi方向,垂线偏差psi方向,水平psi方向,径向
      implicit none
	character*800::loadgrfl
	character*80000::line
	integer nn,n,sn,kln,astat(5)
	real*8::GF(8000,9),rec(8000)
	integer::status=0
!---------------------------------------------------------------------------
      GF=0.d0;nn=0
      open(unit=8,file=loadgrfl,status="old",iostat=status)
      if(status/=0) goto 908
      read(8,'(a)') line
      n=0
      do while(.not.eof(8))  
         read(8,'(a)') line
         call PickReclong(line,kln,rec,sn)
         if(sn<8)goto 906
         rec(2)=rec(2)*1.d-13;rec(3)=rec(3)*1.d-17;rec(4)=rec(4)*1.d-18;rec(5)=rec(5)*1.d-14
         rec(6)=rec(6)*1.d-19;rec(7)=rec(7)*1.d-12;rec(8)=rec(8)*1.d-12;rec(9:10)=rec(9:10)*1.d-21
         n=n+1;GF(n,1:9)=rec(2:10)
906      continue
      enddo
      close(8)
      nn=n
908   continue
	end
!
!****************************************************************************
!
      subroutine IntpGrnF(GF,dl,vfn)
      implicit none
	integer::i,k,kk
	real*8::GF(8000,9),dl,vfn(9),wgh,pr
!---------------------------------------------------------------------------
      kk=nint(dl*1.d-2+0.5d0);if(kk>8000)kk=8000
      pr=0.d0;vfn(1:9)=0.d0
      do i=kk-2,kk+2
        if(i<1.or.i>8000)goto 1100
          wgh=1.d0/(dabs(dble(i)*1.d2-dl)+1.d-10);pr=pr+wgh
          do k=1,9
            vfn(k)=vfn(k)+GF(i,k)*wgh
          enddo
1100    continue
      enddo
      vfn=vfn/pr
      end
