program main

  use iso_fortran_env, only : dp => real64
  use parameters
  use arrays
  use data_types_observables
  use starts
  use matrix_operations
  use dynamics
  implicit none

  integer :: i
  
  call read_input_parameters()
  allocate(U(L,L,L,L))

  call hot_start(U)

  call set_periodic_bounds(L)
  
  !print*, det(U(L,L,L,L)%link(4)%matrix)
  !print*, "U1", U(1,1,1,1)%link(1)
  !print*, "U2", U(1,1,1,1)%link(2)
  !print*, "U1+U2",U(1,1,1,1)%link(1) + U(1,1,1,1)%link(2)
  print*, -action(U,L,1.0_dp/3)/L**4

  do i = 1, 200
     call sweeps(U,L,8.0_dp,3,4,"metropolis")     
     print*, -action(U,L,1.0_dp/3)/L**4
  end do
end program main
