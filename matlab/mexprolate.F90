#include "fintrf.h"

! because prolate subroutines are contained in a module, the %val() construct
! does not work. for now, data copies are used to pass the information into the
! relevant fortran subroutines. to avoid the unnecessary data copies, the code
! will need to be updated to make use of the Matlab APIs in: goo.gl/h2v3I0

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  use prolatemod
  implicit none
  integer, parameter :: dp = SELECTED_REAL_KIND(15, 307)

  ! mex function arguments
  mwPointer :: plhs(*), prhs(*)
  integer :: nlhs, nrhs

  ! function declarations
  mwPointer :: mxGetPr
  mwPointer :: mxCreateDoubleMatrix
  mwPointer :: mxGetM, mxGetN
  integer :: mxIsDouble, mxIsComplex

  ! pointers to input/output mxArrays
  mwPointer :: c_ptr, w_ptr, khi_ptr

  ! array information
  mwPointer :: m, n
  mwSize :: insize, wm, wn, khim, khin

  ! message container for prolate error codes
  character(len=256) msg

  ! arguments for computational routine
  real(dp) :: c
  real(dp), allocatable :: w(:), khi(:)
  integer :: info

  ! check for proper number of arguments
  if ( nrhs /= 1 ) then
     call mexErrMsgIdAndTxt ('MATLAB:prolate:nInput', &
          'One input required.')
  elseif ( nlhs > 2 ) then
     call mexErrMsgIdAndTxt ('MATLAB:prolate:nOutput', &
          'Too many output arguments.')
  endif

  ! validate inputs
  if ( mxIsDouble(prhs(1)) == 0 ) then
     call mexErrMsgIdAndTxt ('MATLAB:prolate:NonDouble', &
          'Input argument must be of type double')
  endif
  if ( mxIsComplex(prhs(1)) == 1 ) then
     call mexErrMsgIdAndTxt ('MATLAB:prolate:NonReal', &
          'Input argument must be a real only')
  endif

  ! get the size of the input array.
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  insize = m*n

  if ( insize /= 1 ) then
     call mexErrMsgIdAndTxt('MATLAB:prolcrea:nInput', &
          'Input must be a scalar')
  end if

  ! create fortran array from the input argument.
  c_ptr = mxGetPr(prhs(1))
  call mxCopyPtrToReal8(c_ptr, c, insize)

  ! call the computational subroutine
  call prolcrea(c, w, khi, info)

  if ( info /= 0 ) then
     write (msg, '(A47,I5)') 'Error in LAPACK subroutine dsteqr, &
          with code = ', info
     call mexErrMsgIdAndTxt('MATLAB:prolate:ProlcreaErr', &
          trim(msg))
  end if

  ! create matrix for the return argument
  ! explicitely set to mwSize variables for compatibility with
  ! -largeArrayDims compilation (necessary for future proofing)
  wm = size(w)
  wn = 1
  khim = size(khi)
  khin = 1

  plhs(1) = mxCreateDoubleMatrix(wm, wn, 0)
  plhs(2) = mxCreateDoubleMatrix(khim, khin, 0)
  w_ptr = mxGetPr(plhs(1))
  khi_ptr = mxGetPr(plhs(2))

  ! load the data into outputs for MATLAB.
  call mxCopyReal8ToPtr(w, w_ptr, wm*wn)
  call mxCopyReal8ToPtr(khi, khi_ptr, khim*khin)
end subroutine mexFunction
