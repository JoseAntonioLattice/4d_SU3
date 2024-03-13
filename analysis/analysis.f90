program analysis

  use iso_fortran_env, only : dp => real64, i4 => int32
  
  implicit none

  integer(i4) :: L
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  character(25) :: algorithm

  namelist /input_parameters/ L, N_thermalization, N_measurements, N_skip, algorithm

  integer(i4) :: inunit
  character(256) :: data_file = "action_metropolis.dat"
  
  type observable
     real(dp), allocatable :: array(:)
     real(dp) :: avr
     real(dp) :: err
  end type observable

  type(observable) :: Ep

  integer(i4) :: i
  
  call read_input()

  N_measurements = 100

  allocate(Ep%array(N_measurements))
  
  data_file = "/home/estudiante/Documents/doctorado/SU3_gauge/data/"//trim(data_file)
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
