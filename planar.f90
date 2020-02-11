program main
  implicit none
  !2D in-plane full-dynamic BIEM with slip weakening friction law.
  !based on Cochard & Madariaga, 1994 and Tada & Madariaga, 2001.
  !developed by SO OZAWA, master student at UTokyo
  real(8),allocatable::kern(:,:),vel(:,:),disp(:),summ(:),tau(:),stress(:)
  integer,allocatable::state(:) !0:locked 1:slipped
  real(8)::dt,dx,t,th,mu,dc,x,tau0,et,tr
  real(8),parameter::cs=1.d0,pi=4.d0*atan(1.d0),cp=cs*sqrt(3d0)
  integer::i,j,k,n,imax,kmax,p,q

  imax=400
  kmax=400
  dx=0.025d0
  dt=dx*0.5d0
  mu=1.d0
  dc=0.1d0
  t=0.d0
  vel=0.d0
  et=0.43d0

  tau0=0.5d0
  th=1.d0
  tr=0.d0

  allocate(kern(-imax+1:imax-1,0:kmax),vel(imax,0:kmax),disp(imax),summ(imax),tau(imax))
  allocate(state(imax),stress(imax))
  kern=0

  open(25,file='kern.dat')

  do p=-imax+1,imax-1
    do q=0,kmax
      kern(p,q)=inte((p+0.5d0)*dx,(q+et)*dt)-inte((p-0.5d0)*dx,(q+et)*dt)-inte((p+0.5d0)*dx,(q-1+et)*dt)+inte((p-0.5d0)*dx,(q-1+et)*dt)
      write(25,'(2i6,4e15.6)') p,q,kern(p,q)
    end do
    write(25,*)
    write(25,*)
  end do
  !stop
  ! do i=1,imax
  !   write(*,*) mu/2.d0/pi*(1.d0/(dx*(i-0.5d0))-1.d0/(dx*(i+0.5d0)))
  ! end do
  ! stop


  call initialcrack(disp,tau,state,imax,dx,th,tau0)
  !do i=1,imax
  !  write(*,*) disp(i),tau(i),state(i)
  !end do
  open(12,file='tmp')
  !disp=0.d0
  ! tau=tau0
  ! do i=250,250
  !   tau(i)=th+0.001d0
  !   state(i)=1
  !   vel(i,0)=0.d0
  ! end do
  !stop

  !time evolution
  !write(*,*) 'start'
  do k=1,kmax
    do i=1,imax
      summ(i)=0.d0
      do j=1,imax
        do n=1,k-1
          summ(i)=summ(i)+vel(j,n)*kern(i-j,k-n)
        end do
      end do
    end do
    !if(vel(i,k)) known, calculate stress directly
    !if(vel(i,k)) unknown, combine with friction law to solve
    do i=1,imax
      if(state(i).eq.1) then
        !write(*,*) k,max(0.d0,th*(1.d0-disp(i)/dc)),tau(i)
      vel(i,k)=cs/cp*(-2.d0*cs/mu*(max(tr,tr+(th-tr)*(1.d0-disp(i)/dc))-tau(i))-summ(i))
    else if(tau(i)-mu/2.d0/cs*summ(i).gt.th) then
      state(i)=1
      vel(i,k)=cs/cp*(-2.d0*cs/mu*(max(tr,tr+(th-tr)*(1.d0-disp(i)/dc))-tau(i))-summ(i))
    else
      vel(i,k)=0.d0
    end if
    if(vel(i,k).lt.0.d0) vel(i,k)=0.d0
    end do

    do i=1,imax
      disp(i)=disp(i)+dt*vel(i,k)
      stress(i)=tau(i)-mu/2.d0/cs*(summ(i)+cp/cs*vel(i,k))
    end do
    !writing output
      !write(12,*) vel(0,k),stress(0),disp(0)
    do i=1,imax
      write(12,'(3i6,4e15.6)') k,i,state(i),vel(i,k),stress(i),disp(i),stress(i)-tau(i)
    end do
    write(12,*)
    !write(12,*)
    t=t+dt
    if(mod(k,10).eq.0) write(*,*) 'time step=',k
    if(maxval(vel(:,k)).le.0.0001d0) then
      write(*,*) k,'slip rate zero'
      exit
    end if
  end do

stop

contains
subroutine kernel(kern,dx,dt,imax,kmax)
  implicit none
  real(8),intent(in)::dx,dt
  integer,intent(in)::imax,kmax
  real(8),intent(inout)::kern(:,:)
  integer::p,q
  open(25,file='kern.dat')
  write(*,*) shape(kern)
  write(*,*) kern(-99,0)
  do p=-imax+1,imax-1
    do q=0,kmax
      write(*,*) p,q,kern(p,q)
      kern(p,q)=inte((p+0.5d0)*dx,(q+et)*dt)-inte((p-0.5d0)*dx,(q+et)*dt)-inte((p+0.5d0)*dx,(q-et)*dt)+inte((p-0.5d0)*dx,(q-et)*dt)
      write(*,*) p,q,kern(p,q)
      write(25,'(2i6,4e15.6)') p,q,kern(p,q)
    end do
    write(25,*)
    write(25,*)
  end do
return
end subroutine kernel

function inte(x,t)
  implicit none
  real(8)::x,t,ss,sp,inte,pa
  real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
  ss=cs*t/x
  sp=cp*t/x
  pa=cs/cp
  inte=0d0
  inte=1d0/pa*theta(x)*theta(t)!+sign(1.d0,x)/pi*theta(t-abs(x)/cp)*&
  !&(-4d0/3d0*pa**3*sqrt(sp**2-1d0)**3+4*pa*(1-pa**2)*sqrt(sp**2-1d0)&
  !&-1d0/pa*acos(abs(x)/cp/t))+sign(1.d0,x)/pi*theta(t-abs(x)/cs)*4d0/3d0&
  !&*sqrt(ss**2-1d0)**3
  if(t-abs(x)/cp.ge.0.d0) then
    inte=inte+sign(1.d0,x)/pi*(-4d0/3d0*pa**3*sqrt(sp**2-1d0)**3+4*pa*(1-pa**2)*sqrt(sp**2-1d0)-acos(abs(x)/cp/t)/pa)
  end if
  if(t-abs(x)/cs.ge.0.d0) then
    inte=inte+sign(1.d0,x)/pi*4d0/3d0*sqrt(ss**2-1d0)**3
  end if
  !write(*,*) inte
  ! if(x.gt.0.d0) then
  !   if((abs(s).le.1.d0)) then
  !     inte=1.d0+1.d0/pi*(sqrt(s**(-2)-1.d0)-acos(s))
  !     !inte=1.d0+1.d0/pi*(sqrt(1.d0-s**2)/s-asin(sqrt(1.d0-s**2))*x/abs(x))
  !   else
  !     inte=1.d0
  !   end if
  ! else
  !   if((abs(s).le.1.d0)) then
  !     inte=1.d0/pi*(sqrt(s**(-2)-1.d0)-acos(s))
  !     !inte=1.d0/pi*(sqrt(1.d0-s**2)/s-asin(sqrt(1.d0-s**2))*x/abs(x))
  !   else
  !   inte=0.d0
  !   end if
  ! end if

  return
end function

real(8) function theta(x)
  implicit none
  real(8):: x
  theta=1.d0
  if(x.lt.0.d0) theta=0.d0
  return
end function

subroutine initialcrack(disp,tau,state,imax,dx,th,tau0)
  implicit none
  integer,intent(in)::imax
  real(8),intent(in)::tau0,th,dx
  real(8),intent(out)::disp(:),tau(:)
  integer,intent(out)::state(:)
  real(8)::lc,xcol,kerns(imax,imax)
  real(8),parameter::pi=4.d0*atan(1.d0)

  lc=dc*4.d0/pi*(th-tr)/(tau0-tr)**2

  disp=0.d0
  state=0
  do i=1,imax
    xcol=dx*(i-imax/2-0.5d0)
    if(xcol**2.lt.lc**2) then
      state(i)=1
      disp(i)=0.67d0*dsqrt(lc**2-xcol**2)
    end if
  end do
  do i=1,imax
    do j=1,imax
      kerns(i,j)=2d0*mu/3.d0/pi*(1.d0/(dx*(i-j-0.5d0))-1.d0/(dx*(i-j+0.5d0)))
    end do
    tau(i)=tau0+sum(disp(:)*kerns(i,:))
  end do

  !gaussian
  state=0
  disp=0d0
  do i=1,imax
    xcol=dx*(i-imax/2-0.5d0)
    if(xcol**2.lt.lc**2) then
      tau(i)=tau0+exp(-xcol**2/lc**2)
    end if
  end do

 !  do i=1,imax
 !    xcol=dx*(i-imax/2-0.5d0)
 !      if(xcol**2.lt.lc**2) then
 !       state(i)=1
 !       tau(i)=1.01d0
 !     end if
 ! end do

  return
end subroutine
end program
