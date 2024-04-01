program main

  use iso_fortran_env, only : dp => real64
  use parameters
  use arrays
  use dynamics
  implicit none

  integer :: i
  
  call read_input_parameters()
  allocate(U(L,L,L,L))
 
  beta = [(i*0.1_dp, i = 1, 80)]
 
  call set_periodic_bounds(L)
  
  call equilibrium_dynamics(U,L,beta,3,4,algorithm,N_thermalization,N_measurements,N_skip,equilibrium)
  
end program main
