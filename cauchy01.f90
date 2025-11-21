program kuramoto_combined
  implicit none

  ! Parameters
  integer, parameter :: N = 500
  integer, parameter :: nstep = 200000          
  integer, parameter :: outer_steps = 201
  real*8, parameter :: pi = 4.d0 * atan(1.d0)
  real*8, parameter :: h = 0.01d0

  ! Variables
  real*8 :: lambda, u, tt, k1, k2
  real*8 :: t, r, s1, aa, bb, theta1, theta2, phi, rea, ima
  real*8 :: dif, diff, pos, neg, sum
  integer :: i, j, l, m, i1, i2

  real*8, allocatable :: omega(:), theta(:), dthetadt(:), yout(:)
  real*8, allocatable :: column_data(:)
  real*8, allocatable :: adj(:,:)
  real*8, allocatable :: gg(:,:), ggg(:,:)
  real*8, allocatable :: AdjTensor(:,:,:)

  complex*16 :: iot, s

  ! For analytical system
  complex*16 :: zp(N), zn(N), W(N), Z(N)
  complex*16 :: k1p(N), k2p(N), k3p(N), k4p(N)
  complex*16 :: k1n(N), k2n(N), k3n(N), k4n(N)
  real*8 :: ord, ord1, count, gamma

  
  integer :: unit_f, unit_b, unit_of, unit_ob

  ! Initialization
  lambda = 1.0d0
  u = 0.0d0
  tt = 0.00d0
  k1 = 1.0d0
  k2 = -0.5d0
  iot = (0.0d0, 1.0d0)
  gamma = 0.05d0
  count = 0.00d0

  call random_seed()
  allocate(omega(N), theta(N), dthetadt(N), yout(N))
  allocate(adj(N,N))
  allocate(column_data(20000))
  allocate(gg(N,100), ggg(N,100))
  allocate(AdjTensor(N,N,outer_steps))

  
  open(newunit=unit_f, file='cauchy_forward01.dat', status='replace', action='write', iostat=i)
  if (i /= 0) stop 'Cannot open cauchy_forward01.dat'
  open(newunit=unit_b, file='cauchy_backward01.dat', status='replace', action='write', iostat=i)
  if (i /= 0) stop 'Cannot open cauchy_backward01.dat'
  open(newunit=unit_of, file='cauchy_ott_forward01.dat', status='replace', action='write', iostat=i)
  if (i /= 0) stop 'Cannot open cauchy_ott_forward01.dat'
  open(newunit=unit_ob, file='cauchy_ott_backward01.dat', status='replace', action='write', iostat=i)
  if (i /= 0) stop 'Cannot open cauchy_ott_backward01.dat'

 !Read natural frequencies from dat file
  open (unit=114, file='cauchy0.05.dat', status='old', access='sequential', action='read')
do m=1, n
   read (114,*) omega(m)
end do


  ! Random initial phases
  do j=1, N
     call random_number(theta1)
     theta(j) = 2.0d0 * pi * theta1
  end do

  ! ============ MAIN LOOP ============
  do i=1, outer_steps
     t = 0.0d0
     pos = 0.0d0
     neg = 0.0d0

     ! Build adjacency matrix
     do i1=1, N
        do i2=1, N
           dif = abs(theta(i2) - theta(i1))
           if (dif >= pi) then
              diff = pi - mod(dif, pi)
           else
              diff = dif
           end if
           if (diff <= (u*pi)) then
              adj(i1,i2) = k1
           else if ((diff >= (u*pi)) .and. (diff <= (tt*pi))) then
              adj(i1,i2) = 1.0d0
           else
              adj(i1,i2) = k2
           end if
        end do
     end do

     AdjTensor(:,:,i) = adj(:,:)

     ! Simulate Kuramoto phases
     do j=1, nstep
        call derivs(lambda, omega, theta, N, dthetadt, adj)
        call rk4(dthetadt, theta, omega, lambda, N, yout, adj)
        do l=1, N
           if (yout(l) < 0.d0) yout(l) = 2.0d0*pi + yout(l)
           theta(l) = mod(yout(l), 2.0d0*pi)
        end do

        aa=0.d0; bb=0.d0; s=(0.d0,0.d0)
        do m=1,N
           theta1 = mod(theta(m), 2.0d0*pi)
           aa = aa + sin(theta1)
           bb = bb + cos(theta1)
           s = s + (adj(1,m) * exp(iot*theta1))
        end do
        s1 = abs(s/(N*1.0d0))
        r = (1.0d0/(N*1.0d0)) * sqrt(aa**2 + bb**2)
        phi = atan2(aa, bb)
        rea = bb/(N*1.0d0)
        ima = aa/(N*1.0d0)
        t = t + h
     end do

     ! Analytical system (zp, zn)
     do l=1,N
        call random_number(theta1)
        call random_number(theta2)
        zp(l) = cmplx(theta1, theta2)
        call random_number(theta1)
        call random_number(theta2)
        zn(l) = cmplx(theta1, theta2)
     end do

     do j=1, nstep
        call compute_W(zp, zn, W, adj)
        call compute_rhs(zp, zn, W, adj, gamma, k1p)
        call compute_rhs(zp + 0.5d0*h*k1p, zn, W, adj, gamma, k2p)
        call compute_rhs(zp + 0.5d0*h*k2p, zn, W, adj, gamma, k3p)
        call compute_rhs(zp + h*k3p, zn, W, adj, gamma, k4p)
        zp = zp + (h/6.d0)*(k1p + 2.d0*k2p + 2.d0*k3p + k4p)

        call compute_rhs(zn, zp, W, adj, gamma, k1n)
        call compute_rhs(zn + 0.5d0*h*k1n, zp, W, adj, gamma, k2n)
        call compute_rhs(zn + 0.5d0*h*k2n, zp, W, adj, gamma, k3n)
        call compute_rhs(zn + h*k3n, zp, W, adj, gamma, k4n)
        zn = zn + (h/6.d0)*(k1n + 2.d0*k2n + 2.d0*k3n + k4n)

        call compute_order_parameter(zp, zn, Z, adj)
     end do

     ord = sum(abs(Z))
     ord1 = sum(abs(W))

     ! Write outputs to the already open files
     if (i <= outer_steps/2) then
        write(unit_f,*) r, s1, tt
        write(unit_of,*) ord/(N*1.0d0), ord1/(N*1.0d0), tt
        if (i < (outer_steps/2 )) then
           u = u + 0.01d0
           tt = tt + 0.01d0
           count = count + 0.01d0
        end if
     else
        write(unit_b,*) r, s1, tt
        write(unit_ob,*) ord/(N*1.0d0), ord1/(N*1.0d0), tt
        if (i > outer_steps/2) then
           u = u - 0.01d0
           tt = tt - 0.01d0
           count = count - 0.01d0
        end if
     end if
  end do

  ! Close files
  close(unit_f)
  close(unit_b)
  close(unit_of)
  close(unit_ob)

  ! Cleanup
  deallocate(omega, theta, dthetadt, yout)
  deallocate(adj, column_data, gg, ggg, AdjTensor)

contains
!Differential equation subroutine
  subroutine derivs(lambda, omega, y, n, dydx, adj)
    implicit none
    integer :: n, i, j
    real*8 :: lambda, aa
    real*8, dimension(:) :: dydx, y, omega
    real*8, dimension(:,:) :: adj
    do i=1,n
       aa = 0.d0
       do j=1,n
          aa = aa + adj(i,j) * sin(y(j)-y(i))
       end do
       dydx(i) = omega(i) + (lambda/(n*1.0d0)) * aa
    end do
  end subroutine derivs

!RK-4 subroutine
  subroutine rk4(dydx, y, omega, lambda, n, yout, adj)
    implicit none
    integer :: n, i
    real*8 :: h, hh, h6, lambda
    real*8, dimension(:) :: dydx, y, yout, omega
    real*8, dimension(:,:) :: adj
    real*8, allocatable :: dym(:), dyt(:), yt(:)

    h = 0.01d0
    hh = h*0.5d0
    h6 = h/6.d0

    allocate(dym(n), dyt(n), yt(n))

    do i=1,n
       yt(i) = y(i) + hh*dydx(i)
    end do
    call derivs(lambda, omega, yt, n, dyt, adj)

    do i=1,n
       yt(i) = y(i) + hh*dyt(i)
    end do
    call derivs(lambda, omega, yt, n, dym, adj)

    do i=1,n
       yt(i) = y(i) + h*dym(i)
       dym(i) = dyt(i) + dym(i)
    end do
    call derivs(lambda, omega, yt, n, dyt, adj)

    do i=1,n
       yout(i) = y(i) + h6*(dydx(i) + dyt(i) + 2.d0*dym(i))
    end do

    deallocate(dym, dyt, yt)
  end subroutine rk4

 !Analytical complex differential equation subroutine 
 subroutine compute_rhs(z1, z2, W, A, gamma, rhs)
    implicit none
    integer :: i
    complex*16, dimension(N) :: z1, z2, W, rhs
    real*8, dimension(N,N) :: A
    real*8 :: gamma
    do i=1,N
       rhs(i) = -gamma*z1(i) + 0.5d0*(conjg(W(i)) - W(i)*(z1(i)**2))
    end do
  end subroutine compute_rhs

!Computing weighted order parameter
  subroutine compute_W(zp, zn, W, A)
    implicit none
    integer :: i, j
    complex*16, dimension(N) :: zp, zn, W
    real*8, dimension(N,N) :: A
    W = (0.d0,0.d0)
    do i=1,N
       do j=1,N
          if (A(i,j) > 0.d0) then
             W(i) = W(i) + A(i,j) * conjg(zp(j))
          else if (A(i,j) < 0.d0) then
             W(i) = W(i) + A(i,j) * conjg(zn(j))
          end if
       end do
       W(i) = W(i)/(N*1.0d0)
    end do
  end subroutine compute_W

!Computing order parameter after stability 
  subroutine compute_order_parameter(zp, zn, Z, A)
    implicit none
    integer :: i, j
    complex*16, dimension(N) :: zp, zn, Z
    real*8, dimension(N,N) :: A
    Z = (0.d0,0.d0)
    do i=1,N
       do j=1,N
          if (A(i,j) > 0.d0) then
             Z(i) = Z(i) + conjg(zp(j))
          else if (A(i,j) < 0.d0) then
             Z(i) = Z(i) + conjg(zn(j))
          end if
       end do
       Z(i) = Z(i)/(N*1.0d0)
    end do
  end subroutine compute_order_parameter

end program kuramoto_combined

