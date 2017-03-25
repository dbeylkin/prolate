program prolate_driver
  use prolatemod
  use printmod
  implicit none
  integer, parameter :: dp = SELECTED_REAL_KIND(15, 307)

  real(dp) :: c
  real(dp), allocatable :: w(:), khi(:)
  integer :: info = 0

  print *, 'Enter c'
  read *, c
  call print('c = ', c)

  call prolcrea(c, w, khi, info)
  if ( info /= 0 ) call print('info = ', info)

  call print('w size = ', size(w))
  call print('khi = ', khi(:))

end program prolate_driver
