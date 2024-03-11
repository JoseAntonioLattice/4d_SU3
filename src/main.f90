program main

  use iso_fortran_env, only : dp => real64
  use parameters
  use arrays
  use data_types_observables
  use starts
  use matrix_operations
  use dynamics
  use statistics
  implicit none

  integer :: i, i_beta
  
  call read_input_parameters()
  allocate(U(L,L,L,L))
  allocate(E_p%array(N_measurements))

  beta = [(i*0.1_dp, i = 1, 80)]
  
  

  call set_periodic_bounds(L)
  
  open(unit = 100, file = 'data/Ep_'//trim(algorithm)//'.dat')

  call hot_start(U)

  do i_beta = 1, size(beta)
     !print*, beta(i_beta), "Before Thermalization"
     call thermalization(U,L,beta(i_beta),3,4,algorithm,N_thermalization)
     !print*, beta(i_beta), "After Thermalization"
     call measurements_sweeps(U,L,beta(i_beta),3,4,algorithm,N_measurements,N_skip,E_p%array)
     !print*, beta(i_beta), "After Measurements"
     call std_err(E_p%array,E_p%avr, E_p%err)
     write(100,*) beta(i_beta),E_p%avr, E_p%err
     print*, beta(i_beta),E_p%avr, E_p%err
  end do
  
  
end program main
