program analysis

  use iso_fortran_env, only : dp => real64, i4 => int32
  use statistics
  use number2string_mod
  implicit none

  integer(i4) :: L
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  character(25) :: algorithm
  logical :: equilibrium
  
  namelist /input_parameters/ L, N_thermalization, N_measurements, N_skip, algorithm, equilibrium

  integer(i4) :: inunit
  character(256) :: directory 
  
  type observable
     real(dp), allocatable :: array(:)
     real(dp) :: avr
     real(dp) :: err
  end type observable

  type(observable) :: Ep

  integer(i4) :: i
  
  call read_input()


  allocate(Ep%array(N_measurements))
  
  data_file = "data/L="//trim(int2str(L))//"/"//trim(eq)//"/"//trim(algorithm)&
       //"/beta="//trim(beta(i_beta))
  print*, trim(data_file)
  open( newunit = inunit, file = trim(data_file) )

  do i = 1, N_measurements
     read(inunit,*) Ep%array(i)
  end do

  print*, Ep%array(1:5)
  
contains

  subroutine read_input()
    integer(i4) :: in_unit
    character(20) :: parameters_file

    read(*,*) parameters_file
    open(newunit = in_unit, file = parameters_file)
    read(in_unit, nml = input_parameters)
    close(in_unit)
    write(*,nml = input_parameters)
    
  end subroutine read_input

  
end program analysis
