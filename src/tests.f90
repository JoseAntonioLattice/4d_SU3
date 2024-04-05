program main
  !use data_types_obersvables
  use matrix_operations
  implicit none

  !integer :: x = 4, d = 3, L = 2

  type(complex_3x3_matrix) :: A,B

  A = generate_SU3()
  B = dagger(A)
  
  print*, det(A%matrix)
  print*, " "
  print*, A*dagger(A)
  print*, " "
  print*, A*B
  print*,matmul(A%matrix,B%matrix)

contains

  function generate_SU3() result(x)
    complex(dp), dimension(3) :: u, v
    type(complex_3x3_matrix) :: x
    real(dp), dimension(6) :: r1, r2
    
    call random_number(r1); call random_number(r2)
    r1 = r1 - 0.5_dp ; r2 = r2 - 0.5_dp;
    r1 = r1/norm2(r1); r2 = r2/norm2(r2)
    u = [cmplx(r1(1),r1(2),dp), cmplx(r1(3),r1(4),dp), cmplx(r1(5), r1(6),dp)]
    v = [cmplx(r2(1),r2(2),dp), cmplx(r2(3),r2(4),dp), cmplx(r2(5), r2(6),dp)]
    v = v - u * (conjg(u(1))*v(1) + conjg(u(2))*v(2) + conjg(u(3))*v(3))!dot_product(u,v)
    v = v / sqrt((v(1)%re)**2 + (v(1)%im)**2 + (v(2)%re)**2 + (v(2)%im)**2 +  (v(3)%re)**2 + (v(3)%im)**2)
    x%matrix = transpose( reshape([u,v,cross(conjg(u), conjg(v))], [3,3]) )
    
  end function generate_SU3
  
  function cross(u,v)
    complex(dp), dimension(3), intent(in) :: u, v
    complex(dp), dimension(3) :: cross

    cross(1) = u(2) * v(3) - u(3) * v(2) 
    cross(2) =-u(1) * v(3) + u(3) * v(1)
    cross(3) = u(1) * v(2) - u(2) * v(1)
  end function cross


  function det(A)

    complex(dp), dimension(3,3), intent(in) :: A
    integer(i4) :: i, info
    integer :: ipiv(3)
    real(dp) :: sgn
    complex(dp) :: det

    ipiv = 0

    call zgetrf(3,3,A,3,ipiv,info)

    det = 1.0_dp

    do i = 1, 3
       det = det*A(i,i)
    end do

    sgn = 1.0_dp
    do i = 1, 3
       if(ipiv(i) /= i) sgn = -sgn
    end do

    det = sgn*det
    
  end function det
  
end program main
