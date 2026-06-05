program exponential
  use iso_fortran_env, only : dp => real64
  implicit none
  integer, parameter :: K = 20
  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: ii = (0.0_dp,1.0_dp)
  !Lie algebra su(3) matrix
  complex(dp), dimension(3,3) :: X, expX, Id
  complex(dp) :: q0, q1, q2, q0old, q1old, q2old
  complex(dp) :: d, t, z 
  
  integer :: i

  print*, pi

  z = 0.5*pi!cmplx(0.0_dp)
  

  Id = 0.0_dp
  Id(1,1) = 1.0_dp
  Id(2,2) = 1.0_dp
  Id(3,3) = 1.0_dp
  
  X = reshape([(0.0_dp,0.744_dp),(-0.192_dp,0.873_dp),(-0.518_dp,0.738_dp), &
       (0.192_dp,0.873_dp),(0.0_dp,-1.444_dp),(0.306_dp,-0.255_dp), &
       (0.518_dp,0.738_dp),(-0.306_dp,-0.255_dp),(0.0_dp,0.7_dp)],shape(X))

  d = ii*determinant(3,X)
  t = -0.5*tr(matmul(X,X))
  print*, d, t
  q0old = 1.0_dp/gamma(1.0_dp*(K+1))
  q1old = (0.0_dp,0.0_dp)
  q2old = (0.0_dp,0.0_dp)
  do i = K-1,0,-1
     q0 = 1.0_dp/gamma(1.0_dp*(i+1)) - ii*d*q2old
     q1 = q0old - t*q2old
     q2 = q1old
     q0old = q0
     q1old = q1
     q2old = q2
     print*, i, q0,q1,q2
  end do

  expX = q0*Id + q1*X + q2*matmul(X,X)

  print*, q0,q1,q2
 ! print*, expX

contains

  function factorial(n)
    integer :: factorial
    integer, intent(in) :: n

    factorial = 1
    do i = 2, n
       factorial = factorial*i 
    end do
  end function factorial

  pure recursive function determinant(n, a) result(det)
    implicit none
    integer, intent(in) :: n
    complex(dp), dimension(n,n), intent(in) :: a
    complex(dp) :: det
    integer :: i, sgn
    complex(dp), dimension(n-1, n-1) :: b
    
    if (n == 1) then
       det = a(1,1)
    else
       det = 0.0
       sgn = 1
       do i = 1, n
          ! Extract submatrix
          b(:, :(i-1)) = a(2:, :i-1)
          b(:, i:) = a(2:, i+1:)
          det = det + sgn * a(1, i) * determinant(n-1, b)
          sgn = -sgn
       end do
    end if
  end function determinant

  function tr(A)
    complex(dp), dimension(:,:), intent(in) :: A
    complex(dp) :: tr
    integer :: i

    tr = 0.0_dp
    do i = 1, size(A(1,:))
       tr = tr + A(i,i)
    end do
  end function tr


end program exponential
