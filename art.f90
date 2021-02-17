program main
  implicit none
  real(8),allocatable::kern(:,:,:),V(:,:),kernf(:,:),pkern(:,:)
  real(8)::dt,dx,t,th,mu,dc,x,tau0,et,tr
  real(8),parameter::cs=1d0,pi=4.d0*atan(1.d0),cp=cs*sqrt(3d0)
  integer::i,j,k,n,imax,kmax,p,q,nm,np1,np2,ns1,ns2,m
  real(8)::plan,r,rn,rp,tmp1,tmp,err
  imax=50
  kmax=100
  dx=0.03d0
  dt=dx*0.5d0/cs
  mu=40.d0
  et=0.5d0

  allocate(kern(imax,imax,0:kmax),V(imax,0:kmax),pkern(10,10),kernf(10,10))
  kern=0d0

  ! do i=1,20
  !   xcol(i)=i*dx
  ! end do
  open(25,file='kern.dat')

  do p=1,50
    do q=1,50
  !p=5
  !q=17
  r=dx*(p-q)
  rn=abs(r)-0.5d0*dx
  rp=abs(r)+0.5d0*dx
  np1=max(0,ceiling(rn/cp/dt)-1)
  np2=max(0,ceiling(rp/cp/dt))
  ns1=max(0,ceiling(rn/cs/dt)-1)
  ns2=max(0,ceiling(rp/cs/dt))
  !write(*,*) np1,np2
  !write(*,*) ns1,ns2
  i=0
      do k=np1,np2+i
        !r=xcol(p)-xcol(q)
        kern(k,p,q)=inte(r+0.5d0*dx,(k+et)*dt)-inte(r-0.5d0*dx,(k+et)*dt)-inte(r+0.5d0*dx,(k-1+et)*dt)+inte(r-0.5d0*dx,(k-1+et)*dt)
      !kern(k,p,q)=inte2(r+0.5d0*dx,(k+et)*dt)-inte2(r-0.5d0*dx,(k+et)*dt)-inte2(r+0.5d0*dx,(k-1+et)*dt)+inte2(r-0.5d0*dx,(k-1+et)*dt)
        !write(25,'(3i6,e15.6)') p,q,k,kern(k,p,q)
      end do
      write(25,'(2i6,e15.6)') p,q,sum(kern(np1:np2+i,p,q))
  !    write(25,*)
  !    write(25,*)
    end do
    write(25,*)
  end do


    k=50
    V=1d0
    do p=1,10
      tmp1=0.d0
      do q=41,50
        r=dx*(p-q)
        rn=abs(r)-0.5d0*dx
        rp=abs(r)+0.5d0*dx
        np1=max(0,ceiling(rn/cp/dt)-1)
        np2=max(0,ceiling(rp/cp/dt))
        !write(*,*)np1,np2
        do m=k-np2,k-np1!k-1
           tmp1=tmp1+V(m,q)*kern(k-m,q,p)
        end do
        kernf(p,q-40)=sum(kern(np1:np2,p,q))
        !end if
      end do
      write(*,*) p,tmp1
    end do

    do k=1,10
      tmp=0d0
      do i=1,10
        if(abs(kernf(i,i))>abs(tmp)) then
          tmp=kernf(i,i)
          imax=i
        end if
      end do
      write(*,*)imax,tmp

      do i=1,10
        do j=1,10
          pkern(i,j)=kernf(i,j)-kernf(i,imax)*kernf(imax,j)/tmp
        end do
        !write(*,'(10e15.5)') pkern(i,:)

      end do
      kernf=pkern
      err=0d0
      do i=1,10
        do j=1,10
          err=err+pkern(i,j)*pkern(i,j)
        end do
      end do
      write(*,*) err
    end do

stop
contains
  function inte(x,t)
    implicit none
    real(8)::x,t,ss,sp,inte,pa
    real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
    ss=cs*t/x
    sp=cp*t/x
    pa=cs/cp
    inte=0d0
    inte=1d0/pa*theta(x)*theta(t)!+sign(1.d0,x)/pi*theta(t-abs(x)/cp)*&
    if(t-abs(x)/cp.ge.0.d0) then
      inte=inte+sign(1.d0,x)/pi*(-4d0/3d0*pa**3*sqrt(sp**2-1d0)**3+4*pa*(1-pa**2)*sqrt(sp**2-1d0)-acos(abs(x)/cp/t)/pa)
    end if
    if(t-abs(x)/cs.ge.0.d0) then
      inte=inte+sign(1.d0,x)/pi*4d0/3d0*sqrt(ss**2-1d0)**3
    end if
    return
  end function

  function inte2(x,t)
    implicit none
    real(8)::x,t,ss,sp,inte2,pa
    real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
    ss=cs*t/x

    inte2=0d0
    if(t-abs(x)/cs.ge.0.d0) then
      inte2=4*pi/abs(x)
    end if
    return
  end function


  real(8) function theta(x)
    implicit none
    real(8):: x
    theta=1.d0
    if(x.lt.0.d0) theta=0.d0
    return
  end function

end program
