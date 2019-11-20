program main
  !$ use omp_lib
  implicit none
  !2D in-plane full-dynamic BIEM with rate and state friction law.
  !developed by SO OZAWA, PhD student at Univ. Tokyo
  include 'mpif.h'

  !elastodynamic convolution kernel for shear and normal stress respectively.
  integer::imax,kmax
  real(8)::dx,dt0,dt,alpha !dt is variable
  real(8),allocatable::kernt(:,:,:),kernn(:,:,:),skernt(:,:),skernn(:,:)
  logical::dynamic

  !physical variables
  real(8),allocatable::V(:),P(:),D(:),S0(:),S(:),N0(:),N(:),rupt(:),Vpast(:,:)
  integer,allocatable::state(:)

  !physical constants
  real(8),parameter::cs=1.d0,pi=4.d0*atan(1.d0),cp=cs*sqrt(3d0),mu=1.d0
  real(8)::ff0,fv0,fp0,fdc,fa,fb,p0,sxx0,sxy0,syy0

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
  integer::i,j,k,m,filesize,rn,nf,q,nstep
  character(128)::filename,command

  dx=0.05d0
  dt=dx*0.5d0
  dt0=dt
  fdc=0.005d0
  et=0.8d0
  alpha=3d-2


  !S0=0.4d0
  !N0=1.d0
  sxx0=1.d0
  sxy0=0.4d0
  syy0=1.d0
  P0=0.5d0
  fv0=1d-9
  fp0=fdc/fv0
  fa=0.016d0
  fb=0.020d0

  !geometry
  !orthogonal
  !xtl(1)=0d0;ytl(1)=0d0;xtr(1)=3d0;ytr(1)=0d0;nc(1)=floor(sqrt((xtl(1)-xtr(1))**2+(ytl(1)-ytr(1))**2)/dx)
  !xtl(2)=2d0;ytl(2)=-2d0;xtr(2)=2d0;ytr(2)=0.5d0;nc(2)=floor(sqrt((xtl(2)-xtr(2))**2+(ytl(2)-ytr(2))**2)/dx)

  !fault geometry by input file
  call get_command_argument(1,filename,status=stat)
  open(30,file=filename)
  read(30,*) nf
  allocate(xtl(nf),xtr(nf),ytl(nf),ytr(nf),nc(nf))
  do i=1,nf
  read(30,*) xtl(i),ytl(i),xtr(i),ytr(i)
  end do
  read(30,*) hypox,hypoy
  !xtl(1)=0d0;ytl(1)=0d0;xtr(1)=7d0;ytr(1)=0d0
  !xtl(2)=3d0;ytl(2)=-0.3d0;xtr(2)=14d0;ytr(2)=-0.3d0
  !xtl(3)=10d0;ytl(3)=0d0;xtr(3)=16d0;ytr(3)=0d0
  !hypox=2.0d0;hypoy=0d0
  !write(*,*) sum(xtl(2:1))
!stop
  do i=1,nf
    nc(i)=floor(sqrt((xtl(i)-xtr(i))**2+(ytl(i)-ytr(i))**2)/dx)
  end do

  !determine the number of elements and time step
  imax=sum(nc)
  kmax=imax*4+1
  if(my_rank.eq.0) write(*,*) 'n,kmax',imax,kmax

  !MPI routine
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

  allocate(kernt(0:kmax,imax,imax_),kernn(0:kmax,imax,imax_),Vpast(0:kmax,imax),V(imax),P(imax),D(imax),summt(imax_),summn(imax_),S0(imax),N0(imax),summtg(imax),summng(imax),rupt(imax))
  allocate(skernt(imax,imax_),skernn(imax,imax_))
  allocate(xcol(imax),ycol(imax),nx(imax),ny(imax),state(imax),S(imax),N(imax))
  allocate(xcol_(imax_),ycol_(imax_),nx_(imax_),ny_(imax_),nm(imax,imax_))
  allocate(xr(imax+1),yr(imax+1),ds(imax),xel(imax),xer(imax),yel(imax),yer(imax))

  !mesh generation
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
      kern11=inte11s(x+0.5d0*ds(j),y)-inte11s(x-0.5d0*ds(j),y)
      kern12=inte12s(x+0.5d0*ds(j),y)-inte12s(x-0.5d0*ds(j),y)
      kern22=inte22s(x+0.5d0*ds(j),y)-inte22s(x-0.5d0*ds(j),y)
      skernt(j,i)=0.5d0*(kern11-kern22)*sin2+kern12*cos2
      skernn(j,i)=-(0.5d0*(kern11+kern22)-0.5d0*(kern11-kern22)*cos2-kern12*sin2)

      do k=0,kmax
        kern11=inte11(x+0.5d0*ds(j),y,(k+et)*dt)-inte11(x-0.5d0*ds(j),y,(k+et)*dt)-inte11(x+0.5d0*ds(j),y,(k-1+et)*dt)+inte11(x-0.5d0*ds(j),y,(k-1+et)*dt)
        kern12=inte12(x+0.5d0*ds(j),y,(k+et)*dt)-inte12(x-0.5d0*ds(j),y,(k+et)*dt)-inte12(x+0.5d0*ds(j),y,(k-1+et)*dt)+inte12(x-0.5d0*ds(j),y,(k-1+et)*dt)
        kern22=inte22(x+0.5d0*ds(j),y,(k+et)*dt)-inte22(x-0.5d0*ds(j),y,(k+et)*dt)-inte22(x+0.5d0*ds(j),y,(k-1+et)*dt)+inte22(x-0.5d0*ds(j),y,(k-1+et)*dt)

        kernt(k,j,i)=0.5d0*(kern11-kern22)*sin2+kern12*cos2
        kernn(k,j,i)=0.5d0*(kern11+kern22)-0.5d0*(kern11-kern22)*cos2-kern12*sin2
        kernt(k,j,i)=kernt(k,j,i)-skernt(j,i)*dt
        kernn(k,j,i)=kernn(k,j,i)-skernn(j,i)*dt

        !if(abs(kernn(i,j,k)).le.1d-8) kernn(k,j,i)=0d0
        !write(25,'(3i6,4e15.6)') i,j,k,kernt(i,j,k),kernn(i,j,k)
      end do
      !write(25,*)
    end do
    !write(25,*)
  end do
  !kernn=-kernn
  !write(*,*)minval(kernt),minval(kernn)

  ! i=50
  ! j=150
  ! k=kmax
  ! if(my_rank.eq.0) then
  !   open(25,file='kern.dat')
  ! do k=1,kmax
  !   write(25,'(3i6,4e15.6)') i,j,k,kernt(k,j,i),kernn(k,j,i),skernt(j,i),skernn(j,i)
  ! end do
  ! end if
  ! stop
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
  t=0d0
  do i=1,imax
    ang=-dasin(nx(i))
    S0(i)=sxy0*cos(2*ang)+0.5d0*(sxx0-syy0)*sin(2*ang)
    N0(i)=sin(ang)**2*sxx0+cos(ang)**2*syy0+sxy0*sin(2*ang)
  end do

  call initialcrack(hypox,hypoy,xcol,ycol,S0,imax,dx)
  P=0.5d0
  Vpast=0d0
  do i=1,imax
    V(i)=fv0*exp((S0(i)/N0(i)-P(i))/fa)
    !write(*,*) i,V(0,i)
  end do
  !stop
  D=0d0
  rupt=0d0
  dynamic=.true.
  k=0
  do nstep=1,6000
    k=k+1
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    time3= MPI_Wtime()
    do i=1,imax_
      summt(i)=0.d0
      summn(i)=0.d0
      !!$omp parallel do
      do j=1,imax
        !if(nm(i,j).lt.k) then
        q=nm(j,i)
        if(dynamic) then
        do m=1,k-1!k-1
          if(Vpast(m,j).ne.0d0) then
           tmp1=Vpast(m,j)*kernt(k-m,j,i)
           tmp2=Vpast(m,j)*kernn(k-m,j,i)
           summt(i)=summt(i)+tmp1
           summn(i)=summn(i)+tmp2
          end if
        end do
        end if
        !static component
        summt(i)=summt(i)+skernt(j,i)*D(j)
        summn(i)=summn(i)+skernn(j,i)*D(j)
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
  !if(my_rank.eq.0) then
    do i=1,imax
      dpdt=fb/fdc*(fv0*exp((p0-P(i))/fb)-V(i))
      P(i)=P(i)+dt*dpdt
      !write(*,*) i,P(k,i)
      !p(k,i)=p(k-1,i)+dt*(1-vel(k-1,i)*p(k-1,i)/dc)
      lnv=rtnewt(dlog(V(i)/fv0),1d-4,N(i),P(i),S0(i),summtg(i))
      Vpast(k,i)=fv0*exp(lnv)
      V(i)=Vpast(k,i)
      if(V(i).gt.0.1d0 .and. state(i).eq.0) then
        state(i)=1
        rupt(i)=t
      end if
    end do

    do i=1,imax
      D(i)=D(i)+dt*V(i)
      S(i)=S0(i)-mu/2.d0/cs*(summtg(i)+cp/cs*V(i))
    end do
    !writing output
      !output
      if(my_rank.eq.0) then
    do i=1,imax
      write(12,'(2i6,7e15.6,i6)') i,k,xcol(i),ycol(i),V(i),S(i),D(i),N(i),S(i)/N(i),state(i)
    end do
    write(12,*)
  end if
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(V, size(V), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    !write(*,*) V(k,:)
    t=t+dt
    time2= MPI_Wtime()
    if(mod(nstep,10).eq.0 .and. my_rank.eq.0) write(*,'(a,i0,3f15.6)') 'time step=',nstep,time2-time1,maxval(V(:)),t
    !if(maxval(V(k,:)).le.1d-4) then
    !  write(*,*) k,'slip rate zero'
    !  exit
    !end if
    !dynamic off
    if(dynamic) then
      dt=dt0
      if(k.gt.kmax) exit
      if(maxval(V(:)).lt.3d-3) dynamic=.false.
      if(k.lt.10) dynamic=.true.
    !dynamic on
    else
      dt=max(dt0,alpha*dt0/maxval(V))
      if(maxval(V(:)).lt.1d-6) exit
      if(maxval(V(:)).gt.1d-2) then
        k=0
        dynamic=.true.
      end if
    end if
  end do

  if(my_rank.eq.0) then
    open(13,file='rupt.dat')
    do i=1,imax
      write(13,*) xcol(i),ycol(i),rupt(i)
    end do
    close(13)
  end if

  !quasi-dynamic simulation

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

function inte12s(x1,x2)
  implicit none
  real(8)::x1,x2,t,ss,sp,inte,pa,inte12s
  real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
  r=sqrt(x1**2+x2**2)
  pa=cs/cp
  inte12s=2*cs/pi*(1-pa**2)*x1*(x1**2-x2**2)/r**4
  return
end function

function inte11s(x1,x2)
  implicit none
  real(8)::x1,x2,t,ss,sp,inte,pa,inte11s
  real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
  r=sqrt(x1**2+x2**2)
  pa=cs/cp
  inte11s=-2*cs/pi*(1-pa**2)*x2*(3*x1**2+x2**2)/r**4
  return
end function

function inte22s(x1,x2)
  implicit none
  real(8)::x1,x2,t,ss,sp,inte,pa,inte22s
  real(8),parameter::pi=4.d0*atan(1.d0),cs=1.d0,cp=sqrt(3.d0)*cs
  r=sqrt(x1**2+x2**2)
  pa=cs/cp
  inte22s=2*cs/pi*(1-pa**2)*x2*(x1**2-x2**2)/r**4
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

function rtnewt(prev,eps,nst,p,t0,sum)
  integer::j
  integer,parameter::jmax=20
  real(8)::rtnewt,prev,eps
  real(8)::f,df,dx,sum,nst,p,t0
  rtnewt=prev
  !write(*,*) rtnewt
  do j=1,jmax
    x=rtnewt
    f=-fv0*exp(x)+cs/cp*(-2*cs/mu*((p+x*fa)*nst-t0)-sum)
    df=-fv0*exp(x)-2*cs**2/cp/mu*fa*nst
    dx=f/df
    rtnewt=rtnewt-dx
    !write(*,*) rtnewt
    if(abs(dx).lt.eps) return
  end do
  write(*,*) 'maximum iteration'
  stop
end function
end program
