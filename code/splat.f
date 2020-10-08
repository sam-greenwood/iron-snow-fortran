
      subroutine splat(y,x,n)
      implicit real*8(a-h,o-z)
      dimension spl(3,12000),x(1),y(1),work(3,12000)
      call drspln(1,n,x,y,spl,work)
      sa=0.0d0
      dum=y(1)
      do 5 i=2,n
      m=i-1
      h=x(i)-x(m)
      sa=(((spl(3,m)*h/4.d0+spl(2,m)/3.d0)*h+spl(1,m)/2.d0)*h+dum)*h+sa
      dum=y(i)
    5 y(i)=sa
      y(1)=0.0d0
      return
      end
