program main

  use iso_fortran_env, only : dp => real64
  use parameters
  use arrays
  use dynamics
  use statistics
  implicit none

  integer :: i
  
  call read_input_parameters()
  allocate(U(L,L,L,L))
 
  beta = [(i*0.1_dp, i = 1, 3)]
 
  call set_periodic_bounds(L)
  
  open(unit = 100, file = 'data/Ep_'//trim(algorithm)//'.dat')

  call equilibrium_dynamics(U,L,beta,3,4,algorithm,N_thermalization,N_measurements,N_skip)
  
end program main
