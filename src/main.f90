program main

  use iso_fortran_env, only : dp => real64
  use parameters
  use arrays
  use dynamics
  implicit none

  integer :: i
  
  call read_input_parameters()
  allocate(U(Lx,Lx,Lx,Lt))
 
  beta = [5.7_dp]![(i*0.1_dp, i = 1, 80)]
 
  call set_periodic_bounds(Lx,Lt)
  
  call equilibrium_dynamics(U,Lx,Lt,beta,3,4,algorithm,N_thermalization,N_measurements,N_skip,equilibrium)
  
end program main
