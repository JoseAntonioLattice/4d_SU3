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

  beta = [(i*0.5, i = 1, 16)]
  
  

  call set_periodic_bounds(L)
  
  !print*, det(U(L,L,L,L)%link(4)%matrix)
  !print*, "U1", U(1,1,1,1)%link(1)
  !print*, "U2", U(1,1,1,1)%link(2)
  !print*, "U1+U2",U(1,1,1,1)%link(1) + U(1,1,1,1)%link(2)
  open(unit = 100, file = 'data/Ep_metropolis.dat')
  !write(100,*) -action(U,L,1.0_dp/3)/L**4


  call hot_start(U)

  do i_beta = 1, size(beta)
  
     do i = 1, N_thermalization
        call sweeps(U,L,beta(i_beta),3,4,"metropolis")     
        !write(100,*) -action(U,L,1.0_dp/3)/L**4
     end do
 
     do i = 1, N_measurements*N_skip
        call sweeps(U,L,beta(i_beta),3,4,"metropolis")     
        if(mod(i,N_skip) == 0)then
           E_p%array(i/N_skip) = -action(U,L,1.0_dp/3)/L**4
           !write(100,*) E_p%array(i/N_skip)
        end if
     end do
     call std_err(E_p%array,E_p%avr, E_p%err)
     write(100,*) beta(i_beta),E_p%avr, E_p%err
  end do
  
  
end program main
