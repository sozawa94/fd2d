program main
  !$ use omp_lib
  implicit none
  !2D in-plane full-dynamic BIEM with rate and state friction law.
  !developed by SO OZAWA, PhD student at Univ. Tokyo
  include 'mpif.h'

  !elastodynamic convolution kernel for shear and normal stress respectively.
  integer::imax,kmax
  real(8)::dt,dx
  real(8),allocatable::kernt(:,:,:),kernn(:,:,:)

  !physical variables
  real(8),allocatable::V(:,:),P(:),D(:),S0(:),S(:),N0(:),N(:),rupt(:)
  integer,allocatable::state(:)

  !physical constants
  real(8),parameter::cs=1.d0,pi=4.d0*atan(1.d0),cp=cs*sqrt(3d0),mu=1.d0
  real(8)::ffp,ffr,fdc
  real(8)::sxx0,syy0,sxy0

  !for fault geometry and nucleation location
  real(8),allocatable::xcol(:),ycol(:),xel(:),xer(:),yel(:),yer(:),nx(:),ny(:),xr(:),yr(:),ds(:)
  real(8),allocatable::xtl(:),xtr(:),ytl(:),ytr(:)
  integer,allocatable::nm(:,:),nc(:)
  real(8),allocatable::xcol_(:),ycol_(:),nx_(:),ny_(:)
  real(8)::hypox,hypoy !hypocenter

  !for MPI communication and performance evaluation
  integer,allocatable::rcounts(:),displs(:)
  integer::ierr,imax_,icomm,np,amari,my_rank,stat
  real(8)::time1,time2,time3,time4

  !others (i.e. temporal memory, variables for loops...)
  real(8),allocatable::summt(:),summn(:),summtg(:),summng(:)
  real(8)::tmp1,tmp2,t,tp,tr,tau0,et,x,y,r,sin2,cos2,kern11,kern12,kern22,up,ur,ang
  real(8)::factor,rad,dpdt,lnv
  integer::i,j,k,m,filesize,rn,nf,q
  character(128)::filename,command

  dx=0.05d0
  dt=dx*0.5d0
  fdc=0.05d0
  et=0.8d0

  ffp=0.7d0
  ffr=0.2d0

  sxx0=1.d0
  sxy0=0.4d0
  syy0=1.d0

  !fault geometry by input file
  call get_command_argument(1,filename,status=stat)
  open(30,file=filename)
  read(30,*) nf
  allocate(xtl(nf),xtr(nf),ytl(nf),ytr(nf),nc(nf))
  do i=1,nf
  read(30,*) xtl(i),ytl(i),xtr(i),ytr(i)
  end do
  read(30,*) hypox,hypoy
!stop
  do i=1,nf
    nc(i)=floor(sqrt((xtl(i)-xtr(i))**2+(ytl(i)-ytr(i))**2)/dx)
  end do


  imax=sum(nc)
  kmax=imax*4+1
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

  allocate(kernt(0:kmax,imax,imax_),kernn(0:kmax,imax,imax_),V(0:kmax,imax),D(imax),summt(imax_),summn(imax_),S0(imax),N0(imax),summtg(imax),summng(imax),rupt(imax))
  allocate(xcol(imax),ycol(imax),nx(imax),ny(imax),state(imax),S(imax),N(imax))
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

  !initial stress as a function of angle
  do i=1,imax
    ang=-dasin(nx(i))
    S0(i)=sxy0*cos(2*ang)+0.5d0*(sxx0-syy0)*sin(2*ang)
    N0(i)=sin(ang)**2*sxx0+cos(ang)**2*syy0+sxy0*sin(2*ang)
  end do

  !nucleation
  call initialcrack(hypox,hypoy,xcol,ycol,S0,imax,dx)
  t=0d0
  V=0d0
  t=0d0
  D=0d0
  rupt=0d0
  time2= MPI_Wtime()
  if(my_rank.eq.0) write(*,*) time2-time1

  do k=1,kmax
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    time3= MPI_Wtime()
    do i=1,imax_
      summt(i)=0.d0
      summn(i)=0.d0
      !!$omp parallel do
      do j=1,imax
        !if(nm(i,j).lt.k) then
        q=nm(j,i)
        !p=1
        do m=1,k-q!k-1
          if(V(m,j).ne.0d0) then
           tmp1=V(m,j)*kernt(k-m,j,i)
           tmp2=V(m,j)*kernn(k-m,j,i)
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
    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(summn,imax_,MPI_REAL8,summng,rcounts,displs,  MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHERv(summt,imax_,MPI_REAL8,summtg,rcounts,displs,  MPI_REAL8,MPI_COMM_WORLD,ierr)

    do i=1,imax
      N(i)=N0(i)+summng(i)
      if(N(i).lt.0.1d0) N(i)=0.1d0
      if(N(i).gt.1.9d0) N(i)=1.9d0
    end do
    !if(vel(i,k)) known, calculate stress directly
    !if(vel(i,k)) unknown, combine with friction law to solve
    do i=1,imax
      tp=ffp*N(i)
      tr=ffr*N(i)
      rad=2*cs**2/cp/mu
      factor=rad*0.5d0*dt*(tp-tr)/fdc
      if(state(i).eq.1) then
        !write(*,*) k,max(0.d0,th*(1.d0-disp(i)/dc)),tau(i)
        if(D(i).le.fdc) then
          V(k,i)=1/(1+factor)*cs/cp*(-2.d0*cs/mu*(tr+(tp-tr)*(1.d0-(D(i)+0.5d0*dt*V(k-1,i))/fdc)-S0(i))-summtg(i))
          V(k,i)=0.5d0*(V(k-1,i)+V(k,i))
        else
          V(k,i)=cs/cp*(-2d0*cs/mu*(tr-S0(i))-summtg(i))
        end if
    else if(S0(i)-mu/2.d0/cs*summtg(i).gt.tp) then
      state(i)=1
      rupt(i)=k*dt
      V(k,i)=cs/cp*(-2.d0*cs/mu*(max(tr,tr+(tp-tr)*(1.d0-D(i)/fdc))-S0(i))-summtg(i))
    else
      V(k,i)=0.d0
    end if
    if(V(k,i).lt.0.d0) V(k,i)=0.d0
    end do

    do i=1,imax
      D(i)=D(i)+dt*V(k,i)
      S(i)=S0(i)-mu/2.d0/cs*(summtg(i)+cp/cs*V(k,i))
    end do
    !writing output
      !output
      if(my_rank.eq.0) then
    do i=1,imax
      write(12,'(2i6,7e15.6,i6)') i,k,xcol(i),ycol(i),V(k,i),S(i),D(i),N(i),S(i)/N(i),state(i)
    end do
    write(12,*)
  end if
    !write(12,*)
    t=t+dt
    time2= MPI_Wtime()
    if(mod(k,10).eq.0 .and. my_rank.eq.0) write(*,*) 'time step=',k,time2-time1,maxval(V(k,:))
    if(maxval(V(k,:)).le.1d-5) then
      write(*,*) k,'slip rate zero'
      exit
    end if
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

subroutine initialcrack(hypox,hypoy,xcol,ycol,S0,imax,dx)
  implicit none
  integer,intent(in)::imax
  real(8),intent(in)::dx,xcol(:),ycol(:),hypox,hypoy
  real(8),intent(out)::S0(:)
  real(8)::lc,kerns(imax,imax),rr,amp
  real(8),parameter::pi=4.d0*atan(1.d0)

  !lc=fdc*4.d0/pi*(up-ur)/(tau0-ur)**2*0.5d0
  lc=1.d0
  !gaussian
   S0=0.5d0
   amp=0.3d0
   do i=1,imax
     !xcol=dx*(i-imax/2-0.5d0)
     rr=(xcol(i)-hypox)**2+(ycol(i)-hypoy)**2
     if(rr.lt.lc**2) then
       S0(i)=S0(i)+amp*exp(-rr/lc**2)
     end if
   end do

  return
end subroutine

! subroutine initialcrack(hypox,hypoy,xcol,ycol,disp,tau,state,imax,dx,up,ur,tau0)
!   implicit none
!   integer,intent(in)::imax
!   real(8),intent(in)::up,ur,tau0,dx,xcol(:),ycol(:),hypox,hypoy
!   real(8),intent(out)::disp(:),tau(:)
!   integer,intent(out)::state(:)
!   real(8)::lc,kerns(imax,imax),rr
!   real(8),parameter::pi=4.d0*atan(1.d0)
!
!   lc=dc*4.d0/pi*(up-ur)/(tau0-ur)**2*0.5d0
!
!  disp=0.d0
!    state=0
!    !hypo=0.5d0
!   do i=1,imax
!     rr=(xcol(i)-hypox)**2+(ycol(i)-hypoy)**2
!     if(rr.lt.lc**2) then
!       state(i)=1
!       disp(i)=0.67d0*dsqrt(lc**2-rr)
!    end if
!   end do
!   do i=1,imax
!     do j=1,imax
!       kerns(i,j)=2d0*mu/3.d0/pi*(1.d0/(dx*(i-j-0.5d0))-1.d0/(dx*(i-j+0.5d0)))
!     end do
!     tau(i)=tau0+sum(disp(:)*kerns(i,:))
!   end do
!
!   !gaussian
!    state=0
!    disp=0d0
!    tau=tau0
!    do i=1,imax
!      !xcol=dx*(i-imax/2-0.5d0)
!      rr=(xcol(i)-hypox)**2+(ycol(i)-hypoy)**2
!      if(rr.lt.lc**2) then
!        tau(i)=up*1.05d0
!        state(i)=1
!      end if
!    end do
!
!  !  do i=1,imax
!  !    xcol=dx*(i-imax/2-0.5d0)
!  !      if(xcol**2.lt.lc**2) then
!  !       state(i)=1
!  !       tau(i)=1.01d0
!  !     end if
!  ! end do
!
!   return
! end subroutine
end program
