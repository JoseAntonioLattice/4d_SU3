module matrix_operations

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables
  implicit none

  complex(dp), dimension(3,3), parameter :: one = reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp],[3,3])
  
  interface operator(+)
    module procedure mat_sum
  end interface

  interface operator(-)
    module procedure mat_sub
  end interface

  interface operator(*)
    module procedure mat_mult
  end interface

  interface operator(/)
     module procedure div_mat_scalar
  end interface operator(/)

contains

    
  pure function mat_sum(a,b) result(c)
    type(complex_3x3_matrix), intent(in) :: a, b
    type(complex_3x3_matrix) :: c
    c%matrix = a%matrix + b%matrix
  end function mat_sum

  pure function mat_sub(a,b) result(c)
    type(complex_3x3_matrix), intent(in) :: a, b
    type(complex_3x3_matrix) :: c
    c%matrix = a%matrix - b%matrix
  end function mat_sub

  pure function mat_mult(a,b) result(c)
    type(complex_3x3_matrix), intent(in) :: a, b
    type(complex_3x3_matrix) :: c
    c%matrix = matmul(a%matrix,b%matrix)
  end function mat_mult

  pure function div_mat_scalar(A,s) result(C)
    type(complex_3x3_matrix), intent(in) :: A
    real(dp), intent(in) :: s
    type(complex_3x3_matrix) :: C

    C%matrix = A%matrix/s
  end function div_mat_scalar
  
  pure function mat_mult2(a,b) result(c)
    type(complex_3x3_matrix), intent(in) :: a, b
    type(complex_3x3_matrix) :: c
    complex(dp), dimension(3) :: cross_prod_b, u, v

    !u = b%matrix(1,:)
    !v = b%matrix(2,:)

    !cross_prod_b = cross_3d(conjg(u),conjg(v))
    
    c%matrix(1,1) = a%matrix(1,1)*b%matrix(1,1) &
                  + a%matrix(1,2)*b%matrix(2,1) &
                  + a%matrix(1,3)*b%matrix(3,1)!cross_prod_b(1)
    c%matrix(1,2) = a%matrix(1,1)*b%matrix(1,2) &
                  + a%matrix(1,2)*b%matrix(2,2) &
                  + a%matrix(1,3)*b%matrix(3,2)
    c%matrix(1,3) = a%matrix(1,1)*b%matrix(1,3) &
                  + a%matrix(1,2)*b%matrix(2,3) &
                  + a%matrix(1,3)*b%matrix(3,3)!cross_prod_b(3)

    c%matrix(2,1) = a%matrix(2,1)*b%matrix(1,1) &
                  + a%matrix(2,2)*b%matrix(2,1) &
                  + a%matrix(2,3)*b%matrix(3,1)!cross_prod_b(1)
    c%matrix(2,2) = a%matrix(2,1)*b%matrix(1,2) &
                  + a%matrix(2,2)*b%matrix(2,2) &
                  + a%matrix(2,3)*b%matrix(3,2)
    c%matrix(2,3) = a%matrix(2,1)*b%matrix(1,3) &
                  + a%matrix(2,2)*b%matrix(2,3) &
                  + a%matrix(2,3)*b%matrix(3,3)!cross_prod_b(3)
    
    c%matrix(3,:) = cross_3d(conjg(c%matrix(1,:)), conjg(c%matrix(2,:)) )

  end function mat_mult2

  pure function cross_3d(u,v) result(w)
    complex(dp), dimension(3), intent(in) :: u, v
    complex(dp), dimension(3) :: w

    w(1) = u(2)*v(3) - u(3)*v(2)
    w(2) = u(3)*v(1) - u(1)*v(3)
    w(3) = u(1)*v(2) - u(2)*v(1)
    
  end function cross_3d
  
  pure function dagger(U) result(U_res)
    type(complex_3x3_matrix), intent(in) :: U
    type(complex_3x3_matrix) :: U_res

    U_res%matrix = transpose(conjg(U%matrix))
  end function dagger

  pure function tr(U) result(trace)
    type(complex_3x3_matrix), intent(in) :: U
    complex(dp) :: trace
    integer :: i

    trace = sum([(U%matrix(i,i),i = 1, size(U%matrix(1,:)))])

  end function tr


end module matrix_operations
