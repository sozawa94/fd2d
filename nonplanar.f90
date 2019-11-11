program main
  !$ use omp_lib
  implicit none
  !2D in-plane full-dynamic BIEM with slip weakening friction law.
  !based on Ando et al., 2007.
  !developed by SO OZAWA, master student at UTokyo
  include 'mpif.h'
  real(8),allocatable::kernt(:,:,:),kernn(:,:,:),vel(:,:),disp(:),summt(:),summn(:),tau(:),stress(:),normal(:),sigma(:),rupt(:)
  real(8),allocatable::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),nx(:),ny(:),xr(:),yr(:),ds(:)
  real(8),allocatable::xcol_(:),ycol_(:),nx_(:),ny_(:),summtg(:),summng(:)
  integer,allocatable::state(:),rcounts(:),displs(:) !0:locked 1:slipping
  integer,allocatable::nm(:,:),nc(:) !1 right 2:left
  real(8)::dt,dx,t,tp,tr,mu,dc,tau0,et,x,y,r,sin2,cos2,kern11,kern12,kern22,up,ur,ang
  real(8),allocatable::xtl(:),xtr(:),ytl(:),ytr(:)
  real(8)::hypox,hypoy,time1,time2,time3,time4,tmp1,tmp2
  real(8),parameter::cs=1.d0,pi=4.d0*atan(1.d0),cp=cs*sqrt(3d0)
  integer::i,j,k,n,imax,kmax,p,q,filesize,rn,nf
  integer::ierr,imax_,icomm,np,amari,my_rank
  character(128)::geofile,command

  dx=0.05d0
  dt=dx*0.5d0
  mu=1.d0
  dc=0.05d0
  et=0.8d0

  tau0=0.4d0
  up=0.6d0
  ur=0.2d0

  !geometry
  !orthogonal
  !xtl(1)=0d0;ytl(1)=0d0;xtr(1)=3d0;ytr(1)=0d0;nc(1)=floor(sqrt((xtl(1)-xtr(1))**2+(ytl(1)-ytr(1))**2)/dx)
  !xtl(2)=2d0;ytl(2)=-2d0;xtr(2)=2d0;ytr(2)=0.5d0;nc(2)=floor(sqrt((xtl(2)-xtr(2))**2+(ytl(2)-ytr(2))**2)/dx)

!step-over
  nf=3
  allocate(xtl(nf),xtr(nf),ytl(nf),ytr(nf),nc(nf))
  xtl(1)=0d0;ytl(1)=0d0;xtr(1)=7d0;ytr(1)=0d0
  xtl(2)=3d0;ytl(2)=-0.3d0;xtr(2)=14d0;ytr(2)=-0.3d0
  xtl(3)=10d0;ytl(3)=0d0;xtr(3)=16d0;ytr(3)=0d0
  hypox=2.0d0;hypoy=0d0
  !write(*,*) sum(xtl(2:1))
!stop
  do i=1,nf
    nc(i)=floor(sqrt((xtl(i)-xtr(i))**2+(ytl(i)-ytr(i))**2)/dx)
  end do


  imax=sum(nc)
  kmax=imax*3+1
  !write(command,"('sed -i -e 's/nstep/',i0,'/g' gnu.gp')") kmax
  !call system(command)
  if(my_rank.eq.0) write(*,*) 'n,kmax',imax,kmax

  !initialize
  icomm=MPI_COMM_WORLD
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr )
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr )
  allocate(rcounts(np),displs(np+1))
  amari=mod(imax,np)
  do k=1,amari
    rcounts(k)=imax/np+1
  end do
  do k=amari+1,np
    rcounts(k)=imax/np
  end do
  displs(1)=0
  do k=2,np+1
    displs(k)=displs(k-1)+rcounts(k-1)
  end do
  imax_=rcounts(my_rank+1)
  write(*,*) imax_

  allocate(kernt(0:kmax,imax,imax_),kernn(0:kmax,imax,imax_),vel(0:kmax,imax),disp(imax),summt(imax_),summn(imax_),tau(imax),normal(imax),summtg(imax),summng(imax),rupt(imax))
  allocate(xcol(imax),ycol(imax),nx(imax),ny(imax),state(imax),stress(imax),sigma(imax))
  allocate(xcol_(imax_),ycol_(imax_),nx_(imax_),ny_(imax_),nm(imax,imax_))
  allocate(xr(imax+1),yr(imax+1),ds(imax),xel(imax),xer(imax),yel(imax),yer(imax))

  do k=1,nf
    do i=0,nc(k)
      xr(i)=xtl(k)+(xtr(k)-xtl(k))*i/nc(k)
      yr(i)=ytl(k)+(ytr(k)-ytl(k))*i/nc(k)
    end do
    do i=1,nc(k)
      j=sum(nc(1:k-1))+i
      xel(j)=xr(i-1)
      xer(j)=xr(i)
      yel(j)=yr(i-1)
      yer(j)=yr(i)
    end do
  end do

  do i=1,imax
    xcol(i)=0.5d0*(xel(i)+xer(i))
    ycol(i)=0.5d0*(yel(i)+yer(i))
    ds(i)=sqrt((xer(i)-xel(i))**2+(yer(i)-yel(i))**2)
    ang=atan((yer(i)-yel(i))/(xer(i)-xel(i)))
    nx(i)=-sin(ang)
    ny(i)=cos(ang)
    !write(*,'(9e15.6)') xcol(i),ycol(i),xel(i),xer(i),yel(i),yer(i),nx(i),ny(i),ds(i)
   end do
   !stop

   call MPI_SCATTERv(xcol,rcounts,displs,MPI_REAL8,xcol_,imax_,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
   call MPI_SCATTERv(ycol,rcounts,displs,MPI_REAL8,ycol_,imax_,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
   call MPI_SCATTERv(nx,rcounts,displs,MPI_REAL8,nx_,imax_,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
   call MPI_SCATTERv(ny,rcounts,displs,MPI_REAL8,ny_,imax_,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  !kernel computation
  time1=MPI_Wtime()
  do i=1,imax_
    do j=1,imax
      x=ny(j)*(xcol_(i)-xcol(j))-nx(j)*(ycol_(i)-ycol(j))
      y=nx(j)*(xcol_(i)-xcol(j))+ny(j)*(ycol_(i)-ycol(j))
      ang=asin(nx_(i))-asin(nx(j))
      sin2=sin(2*ang)
      cos2=cos(2*ang)
      rn=min(sqrt((x+0.5d0*ds(j))**2+y**2),sqrt((x-0.5d0*ds(j))**2+y**2))
      nm(j,i)=max(0,ceiling(rn/cp/dt)-1)
      !write(*,*) nm(i,j)
      !sin2=-2*nx(i)*ny(i)
      !cos2=ny(i)**2-nx(i)**2

      do k=0,kmax
        kern11=inte11(x+0.5d0*ds(j),y,(k+et)*dt)-inte11(x-0.5d0*ds(j),y,(k+et)*dt)-inte11(x+0.5d0*ds(j),y,(k-1+et)*dt)+inte11(x-0.5d0*ds(j),y,(k-1+et)*dt)
        kern12=inte12(x+0.5d0*ds(j),y,(k+et)*dt)-inte12(x-0.5d0*ds(j),y,(k+et)*dt)-inte12(x+0.5d0*ds(j),y,(k-1+et)*dt)+inte12(x-0.5d0*ds(j),y,(k-1+et)*dt)
        kern22=inte22(x+0.5d0*ds(j),y,(k+et)*dt)-inte22(x-0.5d0*ds(j),y,(k+et)*dt)-inte22(x+0.5d0*ds(j),y,(k-1+et)*dt)+inte22(x-0.5d0*ds(j),y,(k-1+et)*dt)

        kernt(k,j,i)=0.5d0*(kern11-kern22)*sin2+kern12*cos2
        kernn(k,j,i)=0.5d0*(kern11+kern22)-0.5d0*(kern11-kern22)*cos2-kern12*sin2
        !if(abs(kernn(i,j,k)).le.1d-8) kernn(k,j,i)=0d0
        !write(25,'(3i6,4e15.6)') i,j,k,kernt(i,j,k),kernn(i,j,k)
      end do
      !write(25,*)
    end do
    !write(25,*)
  end do
  kernn=-kernn
  !write(*,*)minval(kernt),minval(kernn)

  i=50
  j=150
  !k=kmax
  ! if(my_rank.eq.0) then
  !   open(25,file='kern.dat')
  ! do k=1,kmax
  !   write(25,'(3i6,2e15.6)') i,j,k,kernt(k,j,i),kernn(k,j,i)
  ! end do
  ! end if
  !stop
  ! do i=1,imax
  !   write(*,*) mu/2.d0/pi*(1.d0/(dx*(i-0.5d0))-1.d0/(dx*(i+0.5d0)))
  ! end do
  ! stop


  call initialcrack(hypox,hypoy,xcol,ycol,disp,tau,state,imax,dx,up,ur,tau0)
  !do i=1,imax
  !  write(*,*) disp(i),tau(i),state(i)
  !end do
  if(my_rank.eq.0) open(12,file='tmp2')
  !disp=0.d0
  ! tau=tau0
  ! do i=250,250
  !   tau(i)=th+0.001d0
  !   state(i)=1
  !   vel(i,0)=0.d0
  ! end do
  !stop

  !time evolution
  time2= MPI_Wtime()
  if(my_rank.eq.0) write(*,*) time2-time1
  sigma=1d0
  t=0d0
  vel=0d0
  rupt=0d0
  do k=1,kmax
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    time3= MPI_Wtime()
    do i=1,imax_
      summt(i)=0.d0
      summn(i)=0.d0
      !!$omp parallel do
      do j=1,imax
        !if(nm(i,j).lt.k) then
        p=nm(j,i)
        !p=1
        do n=1,k-p!k-1
          if(vel(n,j).ne.0d0) then
           tmp1=vel(n,j)*kernt(k-n,j,i)
           tmp2=vel(n,j)*kernn(k-n,j,i)
           summt(i)=summt(i)+tmp1
           summn(i)=summn(i)+tmp2
          end if
        end do
        !end if
      end do
      !!$omp end parallel do
    end do
    time4= MPI_Wtime()
    !write(*,*) my_rank,time4-time3
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(summn,imax_,MPI_REAL8,summng,rcounts,displs,  MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(summt,imax_,MPI_REAL8,summtg,rcounts,displs,  MPI_REAL8,MPI_COMM_WORLD,ierr)

    do i=1,imax
      normal(i)=sigma(i)+summng(i)
    end do
    !if(vel(i,k)) known, calculate stress directly
    !if(vel(i,k)) unknown, combine with friction law to solve
    do i=1,imax
      tp=up*max(0d0,normal(i))
      tr=ur*max(0d0,normal(i))
      if(state(i).eq.1) then
        !write(*,*) k,max(0.d0,th*(1.d0-disp(i)/dc)),tau(i)
      vel(k,i)=cs/cp*(-2.d0*cs/mu*(max(tr,tr+(tp-tr)*(1.d0-disp(i)/dc))-tau(i))-summtg(i))
    else if(tau(i)-mu/2.d0/cs*summtg(i).gt.tp) then
      state(i)=1
      rupt(i)=k*dt
      vel(k,i)=cs/cp*(-2.d0*cs/mu*(max(tr,tr+(tp-tr)*(1.d0-disp(i)/dc))-tau(i))-summtg(i))
    else
      vel(k,i)=0.d0
    end if
    if(vel(k,i).lt.0.d0) vel(k,i)=0.d0
    end do

    do i=1,imax
      disp(i)=disp(i)+dt*vel(k,i)
      stress(i)=tau(i)-mu/2.d0/cs*(summtg(i)+cp/cs*vel(k,i))
    end do
    !writing output
      !output
      if(my_rank.eq.0) then
    do i=1,imax
      write(12,'(2i6,7e15.6,i6)') i,k,xcol(i),ycol(i),vel(k,i),stress(i),disp(i),normal(i),stress(i)/normal(i),state(i)
    end do
    write(12,*)
  end if
    !write(12,*)
    t=t+dt
    time2= MPI_Wtime()
    if(mod(k,10).eq.0 .and. my_rank.eq.0) write(*,*) 'time step=',k,time2-time1
    ! if(maxval(vel(k,:)).le.0.00001d0) then
    !   write(*,*) k,'slip rate zero'
    !   exit
    ! end if
  end do
  if(my_rank.eq.0) then
    open(13,file='rupt.dat')
    do i=1,imax
      write(13,*) xcol(i),ycol(i),rupt(i)
    end do
    close(13)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

stop

contains
  function inte11(x1,x2,t)
    implicit none
    real(8)::x1,x2,t,ss,sp,inte11,pa,r
    real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
    r=sqrt(x1**2+x2**2)
    ss=cs*t/r
    sp=cp*t/r
    pa=cs/cp
    inte11=0d0
    if(t-abs(r)/cp.ge.0.d0) then
      inte11=inte11-1.d0/pi*2*x2/r*pa*(2*(3*x1**2-x2**2)/(3*r**2)*pa**2*sqrt(sp**2-1d0)**3+(1-2*x2**2/r**2*pa**2)*sqrt(sp**2-1d0))
    end if
    if(t-abs(r)/cs.ge.0.d0) then
      inte11=inte11+1.d0/pi*2*x2/r*(2*(3*x1**2-x2**2)/(3*r**2)*sqrt(ss**2-1d0)**3+(1-2*x2**2/r**2)*sqrt(ss**2-1d0))
    end if
    return
  end function

function inte12(x1,x2,t)
  implicit none
  real(8)::x1,x2,t,ss,sp,inte12,pa,r
  real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
  r=sqrt(x1**2+x2**2)
  ss=cs*t/r
  sp=cp*t/r
  pa=cs/cp
  inte12=0d0
  inte12=1d0/pa*theta(x1)*theta(t-abs(x2)/cs)
  if(t-abs(r)/cp.ge.0.d0) then
    inte12=inte12-sign(1.d0,x1)/pi*2*abs(x1)/r*pa*(2*(3*x2**2-x1**2)/(3*r**2)*pa**2*sqrt(sp**2-1d0)**3+2*x2**2/r*pa**2*sqrt(sp**2-1d0))
  end if
  if(t-abs(r)/cs.ge.0.d0) then
    inte12=inte12+sign(1.d0,x1)/pi*(2*abs(x1)/r*(2*(3*x2**2-x1**2)/(3*r**2)*sqrt(ss**2-1d0)**3+2*x2**2/r**2*sqrt(ss**2-1d0))-acos(abs(x1)/sqrt(cs**2*t**2-x2**2)))
  end if
  return
end function

function inte22(x1,x2,t)
  implicit none
  real(8)::x1,x2,t,ss,sp,inte22,pa,r
  real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
  r=sqrt(x1**2+x2**2)
  ss=cs*t/r
  sp=cp*t/r
  pa=cs/cp
  inte22=0d0
  if(t-abs(r)/cp.ge.0.d0) then
    inte22=inte22+1.d0/pi*2*x2/r*pa*(2*(3*x1**2-x2**2)/(3*r**2)*pa**2*sqrt(sp**2-1d0)**3+(2*x1**2/r**2*pa**2-1d0)*sqrt(sp**2-1d0))
  end if
  if(t-abs(r)/cs.ge.0.d0) then
    inte22=inte22-1.d0/pi*2*x2/r*(2*(3*x1**2-x2**2)/(3*r**2)*sqrt(ss**2-1d0)**3+(2*x1**2/r**2-1d0)*sqrt(ss**2-1d0))
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

subroutine initialcrack(hypox,hypoy,xcol,ycol,disp,tau,state,imax,dx,up,ur,tau0)
  implicit none
  integer,intent(in)::imax
  real(8),intent(in)::up,ur,tau0,dx,xcol(:),ycol(:),hypox,hypoy
  real(8),intent(out)::disp(:),tau(:)
  integer,intent(out)::state(:)
  real(8)::lc,kerns(imax,imax),rr
  real(8),parameter::pi=4.d0*atan(1.d0)

  lc=dc*4.d0/pi*(up-ur)/(tau0-ur)**2*0.5d0

 disp=0.d0
   state=0
   !hypo=0.5d0
  do i=1,imax
    rr=(xcol(i)-hypox)**2+(ycol(i)-hypoy)**2
    if(rr.lt.lc**2) then
      state(i)=1
      disp(i)=0.67d0*dsqrt(lc**2-rr)
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
   tau=tau0
   do i=1,imax
     !xcol=dx*(i-imax/2-0.5d0)
     rr=(xcol(i)-hypox)**2+(ycol(i)-hypoy)**2
     if(rr.lt.lc**2) then
       tau(i)=up*1.05d0
       state(i)=1
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
