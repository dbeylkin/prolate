!*******************************************************************************
!    Fast and accurate algorithm for the computation of prolate spheroidal
!    wave functions (PSWFs).
!    Adapted from prolcrea.f (originally written by V. Rokhlin)
!
!    Copyright (C), 2017 Daniel Beylkin
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*******************************************************************************

module prolatemod
  use printmod
  implicit none
  integer, parameter :: dp = SELECTED_REAL_KIND(15, 307)
  real, parameter :: eps = epsilon(1.0_dp)
  real(dp), parameter :: done = 1
  real(dp), parameter :: half = done/2

  private
  public prolcrea

contains

  subroutine prolcrea(c, w, khi, info)
    real(dp), intent(in) :: c
    real(dp), allocatable, intent(out) :: w(:), khi(:)
    integer, intent(out) :: info

    integer :: nvectss(20), ns(20)
    data nvectss/ 25, 34, 43, 50, 58, 65, 72, 80, 87, 94, 100, &
                  107, 114, 121, 128, 134, 141, 148, 154, 161 /
    data ns/ 48, 64, 80, 92, 106, 120, 130, 144, 156, 168, &
             178, 190, 202, 214, 224, 236, 248, 258, 268, 280 /

    ! local variables
    integer :: i, n, nvects, nhigh, iscale

    ! scratch memory
    real(dp), allocatable :: work(:)

    ! work array memory management variables
    integer :: iladdrw, lladdrw, istorew, lstorew, ltotw
    integer :: istart, iladdr, lladdr, istore, lstore, ltot

    ! initialize status
    info = 0

    n = int(c*3/2)
    nvects = int(c)

    i = int(c/10)
    if ( i <= 19 ) n = ns(i+1)
    if ( i <= 19 ) nvects = nvectss(i+1)

    n = n + 30
    nvects = nvects + 29
    call print('n = ', n)
    call print('nvects = ', nvects)

    ! allocate memory for the subroutine prolvect that obtains coefficients
    ! of Legendre expansions of prolate spheroidal wave functions

    ! allocate memory for work array
    iladdrw = 1
    lladdrw = nvects*4 + 6

    istorew = iladdrw + lladdrw
    lstorew = n*n/2 + n + 6

    ltotw = istorew + lstorew

    if ( allocated(khi) ) deallocate(khi)
    allocate( khi(nvects) )

    allocate ( work(ltotw) )

    call prolvect(n, c, work(istorew), work(iladdrw), khi, nhigh, iscale, info)
    if ( info /= 0 ) return

    ! allocate memory for output array w
    istart = 101

    iladdr = istart
    lladdr = nvects*4

    istore = iladdr + lladdr
    lstore = istore + iscale + nhigh - 2

    ltot = istore + lstore

    if ( allocated(w) ) deallocate(w)
    allocate( w(ltot) )

    ! copy necessary data over to array w
    do i = 1, lladdr
       w(iladdr+i-1) = work(iladdrw+i-1)
    end do

    do i = 1, lstore
       w(istore+i-1) = work(istorew+i-1)
    end do

    ! store varous types of integer data in the beginning of array w
    w(1) = iladdr
    w(2) = istore
    w(3) = istore + iscale - 1
    w(4) = nhigh

    call print('nhigh = ', nhigh)

  end subroutine prolcrea





  subroutine prolvect(n, c, store, laddr, rlamouts, nhigh, iscale, info)
    integer, intent(in) :: n
    real(dp), intent(in) :: c
    real(dp), intent(out) :: store(*), laddr(4,*)
    real(dp), intent(out) :: rlamouts(:)
    integer, intent(out) :: nhigh, iscale, info

    ! scratch memory for eigenvalue solver routines
    real(dp), allocatable :: as(:), bs(:), cs(:), &
         zs1(:,:), zs2(:,:), work(:)

    integer :: nvects, i, j, istore, ihigh
    integer :: infop = 0

    nvects = size(rlamouts)

    ! over allocate slightly for prolmatr
    allocate ( as(n/2+6) )
    allocate ( bs(n/2+6) )
    allocate ( cs(n/2+6) )

    allocate ( zs1(n/2, n/2) )
    allocate ( zs2(n/2, n/2) )
    allocate ( work(n) )

    ! construct the tridiagonal matrix whose eigenvalues are the
    ! "separation coefficients" for the prolate spheroidal wave functions

    ! for the even-numbered separation coefficients
    call prolmatr(as, bs, cs, n, c, 0._dp, .true., .false.)

    ! find the spectrum of the tridiagonal matrix
    ! note db: previous version of the code used a symmetric tridiagonal matrix
    ! solver where the off-diagonal array has entries 2:n/2, i.e., e(1) is
    ! arbitrary. lapack solvers expects array entries 1:n/2-1, i.e., e(n/2)
    ! is arbitrary. since memory is slightly over-allocated, make the adjustment
    ! here since it is simpler than adjusting prolmatr
    call dsteqr('I', n/2, bs, as(2), zs1, n/2, work, infop)
    if ( infop /= 0 ) info = infop
    if ( info /= 0 ) return

    j = 1
    do i = 1, n/2
       rlamouts(j) = -bs(n/2-i+1)
       j = j+2
       if ( j > nvects ) exit
    end do

    ! for the odd-numbered separation coefficients
    call prolmatr(as, bs, cs, n, c, 0._dp, .true., .true.)

    ! find the spectrum of the tridiagonal matrix
    call dsteqr('I', n/2, bs, as(2), zs2, n/2, work, infop)
    if ( infop /= 0 ) info = infop
    if ( info /= 0 ) return

    j = 2
    do i = 1, n/2
       rlamouts(j) = -bs(n/2-i+1)
       j = j+2
       if ( j > nvects ) exit
    end do

    ! store non-zero parts of zs1, zs2 in the array store
    istore = 1
    j = 1
    do i = 1, n/2
       call prolpack(j, zs1(:,n/2-i+1), store, laddr(:,j), istore)
       j = j+1
       if ( j > nvects) exit
       call prolpack(j, zs2(:,n/2-i+1), store, laddr(:,j), istore)
       j = j+1
       if ( j > nvects) exit
    end do

    ! find the highest order of the Legendre expansion for any eigenfunction
    nhigh = 1
    do  i = 1,nvects
       ihigh = int(laddr(3,i)) + int(laddr(4,i))
       if( ihigh > nhigh) nhigh = ihigh
    end do
    nhigh = 2*nhigh + 2

    ! construct and put in arrays store the scaling coefficients
    ! to be used by the subroutine prolev0
    iscale = istore
    do i = 1, nhigh
       store(iscale+i-1) = sqrt(i-half)
    end do

  end subroutine prolvect





  subroutine prolmatr(as, bs, cs, n, c, rlam, ifsymm, ifodd)
    real(dp), intent(in) :: c, rlam
    integer, intent(in) :: n
    logical, intent(in) :: ifsymm, ifodd
    real(dp), intent(out) :: as(*), bs(*), cs(*)

    integer :: k, k0
    real(dp) :: alpha0, beta0, gamma0

    if ( ifodd ) then
       ! construct the tridiagonal matrix corresponding
       ! to odd-numbered P_k
       k = 0
       do k0 = 1, n+2, 2
          k = k+1
          call prolcoef(rlam, k0, c, alpha0, beta0, gamma0, &
               as(k), bs(k), cs(k))

          ! remembering that the norm of P_n is not equal to 1,
          ! rescale the matrix to make it symmetric
          if ( .not. ifsymm ) cycle

          if ( k0 > 1 ) as(k) = as(k)/sqrt(k0-2+half)*sqrt(k0+half)
          cs(k) = cs(k)*sqrt(k0+half)/sqrt(k0+half+2)
       end do
    else
       ! construct the tridiagonal matrix corresponding
       ! to even-numbered P_k
       k = 0
       do k0 = 0, n+2, 2
          k = k+1
          call prolcoef(rlam, k0, c, alpha0, beta0, gamma0, &
               as(k), bs(k), cs(k))

          ! remembering that the norm of P_n is not equal to 1,
          ! rescale the matrix to make it symmetric
          if ( .not. ifsymm ) cycle

          if ( k0 /= 0 ) as(k) = as(k)/sqrt(k0-2+half)*sqrt(k0+half)
          cs(k) = cs(k)*sqrt(k0+half)/sqrt(k0+half+2)
       end do
    end if

  end subroutine prolmatr





  !*****************************************************************************
  ! prolcoef evaluates the legendre coefficients
  ! alpha0, beta0, gamma0, alpha, beta, gamma of two functions:
  !
  !       (1-x**2)   P_k (x)  =
  !                                                                     (1)
  !       alpha0 * P_{k-2} + beta0 * P_{k} + gamma0 * P_{k+2},
  !
  !    and
  !
  !       D               D
  !       -  ( (1-x**2) * -  P_k (x) ) + (rlam-c**2*x**2) * P_k(x) =
  !       Dx              Dx
  !                                                                     (2)
  !       alpha * P_{k-2} + beta * P_{k} + gamma * P_{k+2}.
  !
  !                      input parameters:
  !
  !  rlam - the coefficient (real) in the formula (2)
  !  k - the index in the formulae (1), (2)
  !  c - the coefficient (real) in the formula (2)
  !
  !                      output parameters:
  !
  !  alpha0, beta0, gamma0 - coefficients in the expansion (1)
  !  alpha, beta, gamma - coefficients in the expansion (2)
  !*****************************************************************************
  subroutine prolcoef(rlam, k, c, alpha0, beta0, gamma0, &
                      alpha, beta, gamma)
    integer, intent(in) :: k
    real(dp), intent(in) :: rlam, c
    real(dp), intent(out) :: alpha0, beta0, gamma0, alpha, beta, gamma
    real(dp) :: d, d2, uk, vk, wk

    d = k*(k-1)
    d = d/(2*k+1)/(2*k-1)
    uk = d

    d = (k+1)**2
    d = d/(2*k+3)
    d2 = k**2
    d2 = d2/(2*k-1)
    vk = (d+d2)/(2*k+1)

    d = (k+1)*(k+2)
    d = d/(2*k+1)/(2*k+3)
    wk = d

    alpha = -c**2*uk
    beta = rlam-k*(k+1)-c**2*vk
    gamma = -c**2*wk

    alpha0 = uk
    beta0 = vk
    gamma0 = wk

  end subroutine prolcoef





  subroutine prolpack(k, z, store, laddr, istore)
    integer, intent(in) :: k
    real(dp), intent(in) :: z(:)
    real(dp), intent(out) :: store(*), laddr(4,1)
    integer, intent(inout) :: istore

    integer :: n, i, i1, i2, nn
    n = size(z)

    ! find the first element of the vector z that is non-zero
    i1 = 1
    do i = 1, n
       i1 = i
       if ( abs(z(i)) < eps ) cycle
       exit
    end do

    ! find the last coefficient of the vector z that is non-zero
    i2 = n
    do i = n, i1, -1
       i2 = i
       if ( abs(z(i)) < eps ) cycle
       exit
    end do

    ! store the non-zero elements of this vector (and negate them)
    nn = i2 - i1 + 1
    do i = 1, nn
       store(istore+i-1) = -z(i1+i-1)
    end do

    ! enter the appropriate information in the array laddr
    laddr(1, 1) = k
    laddr(2, 1) = istore
    laddr(3, 1) = nn
    laddr(4, 1) = i1
    istore = istore + nn
  end subroutine prolpack

end module prolatemod

