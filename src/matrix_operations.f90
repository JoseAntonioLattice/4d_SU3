module matrix_operations

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables
  implicit none

  interface operator(+)
    module procedure mat_sum
  end interface

  interface operator(-)
    module procedure mat_sub
  end interface

  interface operator(*)
    module procedure mat_mult
  end interface


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
