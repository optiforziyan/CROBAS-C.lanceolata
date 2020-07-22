!********************************************CROBAS Models for C.lanceolata-Ziyan Liao(2018)****************************************************!
program Crobas
 implicit none 
 real*8 hs,wf,n0,f(25),m(29),i(28),k(28),h,v,hc,d,nfinal,C,hb,b(44)
  !Input data:D,H,N;Calculate data:wf,hs,hc,hb,v,c
  !array 300 means minimum plot numbers.
  !f(25),m(29),i(28),k(28)are symbols;b(44) represent 44 parameters.;
 integer*4 ii,J,th
 real*8, parameter ::  pi=3.1415926
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Set Parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
 write(*,*)'CROBAS Models for C.lanceolata-2018'
  b(2)=1.0_8 !unchange
  b(3)=0.75_8 !unchange
  B(4)=1_8 ! 0-1 ,As this value increases,the DBH decreses.
  b(5)=1.0_8 !unchange
  b(6)=0.23_8 !unchange
  b(7)=1_8 !unchange
  b(8)=300_8 !unchange
  b(9)=0.001_8    !very important,0-0.01,As this value increases,the DBH decreses.
  b(10)=0.00068_8 !
  b(11)=0.00035_8 !
  b(12)=0.0175_8  ! 0.0175~0.2399
  b(13)=2.499_8 !
  b(14)=0.0492_8 !
  b(15)=0.75_8  !0.55~0.95
  b(16)=0.2_8 !unchange
  b(17)=0.02_8 !unchange
  b(18)=0.25_8 !unchange
  b(19)=1.0_8 !unchange
  b(20)=1.0_8 !unchange
  b(21)=0.01_8 !guess
  b(22)=1.0_8 !unchange
  b(23)=0.5_8  !guess
  b(24)=0.9_8 !guess
  b(25)=0.0_8 !unchange
  b(26)=4_8 ! 5.34~6.17
  b(27)=2.134_8 ! important 2.4~9.2
  b(28)=0.02_8 !unchangge
  b(29)=0.2_8  !0.2~0.5
  b(30)=1.6_8 !try 1-10
  b(31)=9.7914_8 !try 1-10
  b(32)=0.03912_8 !try 0-1 As this value increases,the DBH decreses.
  b(33)=0.001_8 !no use for planted forest
  b(34)=0.01_8 !no use for planted forest
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Input Parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !write(*,*)"Please Set The Initial value of X(1):"
  !read(*,*) b(15)
  !write(*,*)"Please Set The Initial value of X(2):"
  !read(*,*) b(26)
  !write(*,*)"Please Set The Initial value of X(3):"
  !read(*,*) b(27)
  !write(*,*)"Please Set The Initial value of X(4):"
  !read(*,*) b(29)
  !write(*,*)"Please Set The Initial value of X(5):"
  !read(*,*) b(32)
  !write(*,*)"Please Set The Initial value of X(6):"
  !read(*,*) b(33)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Other Parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  b(35)=b(6)*b(4)
  b(36)=b(7)*b(5)
  b(37)=b(24)*b(6)
  b(38)=b(25)*b(7)
  b(39)=b(8)*b(9)*b(2)+b(8)*b(11)*b(36)
  b(40)=b(8)*b(9)*b(3)+b(8)*b(10)*b(35)+b(8)*b(11)*b(36)
  b(41)=b(8)*b(9)*(b(2)+b(20)*b(23))+b(8)*b(10)*b(37)*b(20)+b(8)*b(11)*(b(36)+b(38)*b(20))
  b(42)=b(8)*b(9)*b(22)*b(20)+b(8)*b(11)*b(38)*b(20)
  b(43)=b(8)*b(9)*b(22)*b(21)+b(8)*b(11)*b(38)*b(21)
  b(44)=b(8)*b(9)*b(23)*b(21)+b(8)*b(10)*b(37)*b(21)+b(8)*b(11)*b(25)*b(21)
  f(8)=0.0
  m(16)=0.0
  m(29)=0.0
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Record datas~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
 open(unit=111,file='results.csv',status='unknown')
 write(111,*)  'wf',',','hs',',','n0',',','h',',','hc',',','hb',',','v',',','d',',','c'
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Input plots data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
 write(*,*)"Please Set The Initial value of DBH:"
  read(*,*) d
 write(*,*)"Please Set The Initial value of Stand Height:"
  read(*,*) h
 write(*,*)"Please Set The Initial value of Stand Density:"
  read(*,*) n0
 wf= 0.7906*dexp(0.1322*d) !3.3743d0
 hs= 2.347*dexp(0.0743*h)      !5.597351d0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Predicting how many years ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
 write(*,*)"Please Set The Initial value of predicting year:"
  read(*,*) j
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~models~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
 if (n0>1500) then
 do ii=0,j
 hs=hs+f(8)
 wf=wf+m(16)
 n0=abs(n0+m(29))
 !state of tree
 i(2)=wf
 i(3)=b(12)*wf  !wr
 i(8)=hs
 i(7)=(wf/b(14))**(1.0/b(13))  !hc
 i(4)=b(8)*b(9)*(b(2)*hs+b(3)*i(7))*wf  !ws
 i(5)=b(8)*b(10)*b(4)*b(6)*i(7)*wf   !wb
 i(6)=b(8)*b(11)*b(5)*b(7)*(hs+i(7))*wf  !wt
 i(9)=b(6)*i(7)  !hb
 i(10)=b(7)*(i(7)+hs)  !ht
 ! state of stand
 k(3)=abs(n0)/10000.0 !N
 k(2)=b(26)*k(3)*wf !LAI
 k(4)=k(3)*pi*i(9)**2 !C
 if ((b(32)*k(4))**b(30)<1.0)  then
 f(14)=(b(32)*k(4))**b(30) !s
 else
 f(14)=1.0
 end if
 f(16)=b(27)*(1.0d0/k(3))*(1.0d0-exp(-b(29)*k(2))) ! P 光合作用
 f(15)=f(16)/wf  !delta c0
 f(11)=f(15)*(1.0d0-b(28)*i(7))  !delta c
 i(27)=f(11)*wf  
 i(26)=b(16)*(i(3)+wf)+b(17)*(i(4)+i(5)+i(6))  !Rm
 m(28)=((1.0d0)/b(15))*(f(11)-b(16)*(1+b(12))-b(17)*(b(39)*i(8)+b(40)*i(7)))  !Pn
 m(26)=(m(28)-b(18)-b(12)*b(19)-b(43)*hs-b(44)*i(7))/((1.0d0+b(12))+b(39)*hs+((b(13)+1.0)/b(13))*b(40)*i(7))  !R0
 m(27)=(1.0d0-f(14))*m(26)  ! rf(us=0)
 f(17)=-(b(33)+b(34)*k(4)**b(31))*k(3)
 m(29)=f(17)*10000.0d0
 !growth
 k(16)=b(18)*wf  !sf
 k(17)=b(19)*i(3)
 f(8)=(m(28)-b(18)-b(12)*b(19)-b(43)*hs-b(44)*i(7)-(m(27)*(((1.0d0+b(12))+b(39)*hs+((b(13)+1.0d0)/b(13))*b(40)*i(7)))))/ &
  (b(41)+b(42)*(hs/i(7)))  !us
 f(7)=(m(28)-b(18)-b(12)*b(19)-(b(41)+b(42)*(hs/i(7)))*f(8)-b(43)*hs-b(44)*i(7))/((1.0d0+b(12))+b(39)*hs+((b(13)+1.0d0)/ &
  b(13))*b(40)*i(7))
 k(26)=b(20)*(f(8)/i(7))+b(21)
 k(27)=b(20)*(f(8)/i(7))+b(21)
 k(28)=b(20)*(f(8)/i(7))+b(21)
 k(18)=b(8)*b(9)*k(26)*wf*(b(22)*hs+b(23)*i(7))
 k(19)=b(8)*b(10)*b(37)*wf*(b(20)*f(8)+b(21)*i(7))
 k(20)=b(8)*b(11)*k(28)*wf*b(38)*(hs+i(7))
 k(21)=k(16)+k(17)+k(18)+k(19)+k(20)
 i(16)=(f(7)+b(18))*wf
 i(17)=b(12)*(f(7)+b(19))*wf
 m(16)=f(7)*wf
 m(17)=b(12)*f(7)*wf
 m(18)=b(8)*b(9)*(b(2)*f(8)+b(2)*hs*f(7)+f(7)*b(3)*((b(13)+1.0d0)/b(13))*i(7))*wf
 m(19)=b(8)*b(10)*wf*(f(7)*b(35)*((b(13)+1.0d0)/b(13))*i(7))
 m(20)=b(8)*b(11)*b(36)*(f(8)+f(7)*hs+f(7)*((b(13)+1.0d0)/b(13))*i(7))*wf
 m(21)=m(16)+m(17)+m(18)+m(19)+m(20)
 i(18)=k(18)+m(18)  
 i(19)=k(19)+m(19)
 i(20)=k(20)+m(20)
 i(21)=i(16)+i(17)+i(18)+I(19)+i(20)
 i(28)=i(27)-(i(21)/(1.0d0/b(15)))
 hc=i(7)
 h=hs+hc
 C=k(4)
 d=log(wf/0.7906)/0.1322
 v=((h+3)*d**2*pi*0.42*n0)/40000.0
 hb=i(9)
 write(111,'(f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') wf,',',hs,',',n0,',',h,',',hc,',',hb,',',v,',',d,',',c
 write(*,"(T3,i3,'-year')") ii
 write(*,"(f10.2)") h,d
 end do
 else
 do ii=0,j
 hs=hs+f(8)
 wf=wf+m(16)
 n0=n0
 !state of tree
 i(2)=wf
 i(3)=b(12)*wf  !wr
 i(8)=hs
 i(7)=(wf/b(14))**(1.0/b(13))  !hc
 i(4)=b(8)*b(9)*(b(2)*hs+b(3)*i(7))*wf  !ws
 i(5)=b(8)*b(10)*b(4)*b(6)*i(7)*wf   !wb
 i(6)=b(8)*b(11)*b(5)*b(7)*(hs+i(7))*wf  !wt
 i(9)=b(6)*i(7)  !hb
 i(10)=b(7)*(i(7)+hs)  !ht
 ! state of stand
 k(3)=abs(n0)/10000.0 !N
 k(2)=b(26)*k(3)*wf !LAI
 k(4)=k(3)*pi*i(9)**2 !C
 if ((b(32)*k(4))**b(30)<1.0)  then
 f(14)=(b(32)*k(4))**b(30) !s
 else
 f(14)=1.0
 end if
 f(16)=b(27)*(1.0d0/k(3))*(1.0d0-exp(-b(29)*k(2))) ! P 光合作用
 f(15)=f(16)/wf  !delta c0
 f(11)=f(15)*(1.0d0-b(28)*i(7))  !delta c
 i(27)=f(11)*wf  
 i(26)=b(16)*(i(3)+wf)+b(17)*(i(4)+i(5)+i(6))  !Rm
 m(28)=((1.0d0)/b(15))*(f(11)-b(16)*(1+b(12))-b(17)*(b(39)*i(8)+b(40)*i(7)))  !Pn
 m(26)=(m(28)-b(18)-b(12)*b(19)-b(43)*hs-b(44)*i(7))/((1.0d0+b(12))+b(39)*hs+((b(13)+1.0)/b(13))*b(40)*i(7))  !R0
 m(27)=(1.0d0-f(14))*m(26)  ! rf(us=0)
 f(17)=-(b(33)+b(34)*k(4)**b(31))*k(3)
 m(29)=f(17)*10000.0d0
 !growth
 k(16)=b(18)*wf  !sf
 k(17)=b(19)*i(3)
 f(8)=(m(28)-b(18)-b(12)*b(19)-b(43)*hs-b(44)*i(7)-(m(27)*(((1.0d0+b(12))+b(39)*hs+((b(13)+1.0d0)/b(13))*b(40)*i(7)))))/ &
  (b(41)+b(42)*(hs/i(7)))  !us
 f(7)=(m(28)-b(18)-b(12)*b(19)-(b(41)+b(42)*(hs/i(7)))*f(8)-b(43)*hs-b(44)*i(7))/((1.0d0+b(12))+b(39)*hs+((b(13)+1.0d0)/ &
  b(13))*b(40)*i(7))
 k(26)=b(20)*(f(8)/i(7))+b(21)
 k(27)=b(20)*(f(8)/i(7))+b(21)
 k(28)=b(20)*(f(8)/i(7))+b(21)
 k(18)=b(8)*b(9)*k(26)*wf*(b(22)*hs+b(23)*i(7))
 k(19)=b(8)*b(10)*b(37)*wf*(b(20)*f(8)+b(21)*i(7))
 k(20)=b(8)*b(11)*k(28)*wf*B(38)*(hs+i(7))
 k(21)=k(16)+k(17)+k(18)+k(19)+k(20)
 i(16)=(f(7)+b(18))*wf
 i(17)=b(12)*(f(7)+b(19))*wf
 m(16)=f(7)*wf
 m(17)=b(12)*f(7)*wf
 m(18)=b(8)*b(9)*(b(2)*f(8)+b(2)*hs*f(7)+f(7)*b(3)*((b(13)+1.0d0)/b(13))*i(7))*wf
 m(19)=b(8)*b(10)*wf*(f(7)*b(35)*((b(13)+1.0d0)/b(13))*i(7))
 m(20)=b(8)*b(11)*b(36)*(f(8)+f(7)*hs+f(7)*((b(13)+1.0d0)/b(13))*i(7))*wf
 m(21)=m(16)+m(17)+m(18)+m(19)+m(20)
 i(18)=k(18)+m(18)  
 i(19)=k(19)+m(19)
 i(20)=k(20)+m(20)
 i(21)=i(16)+i(17)+i(18)+I(19)+i(20)
 i(28)=i(27)-(i(21)/(1.0d0/b(15)))
 hc=i(7)
 C=k(4)
 d=log(wf/0.7906)/0.1322
 h=0.937*d-3.6298
 v=((h+3)*d**2*pi*0.42*n0)/40000.0
 hb=i(9)
 write(111,'(f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2,a,f10.2)') wf,',',hs,',',n0,',',h,',',hc,',',hb,',',v,',',d,',',c
 write(*,"(T3,i3,'-year')") ii
 write(*,"(f10.2)") h,d
 end do
 end if
end program Crobas
 
!********************************************CROBAS Models for C.lanceolata-Ziyan Liao(2018)****************************************************!