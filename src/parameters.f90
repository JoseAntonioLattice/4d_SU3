module parameters

  use iso_fortran_env, only: dp => real64, i4 => int32

  implicit none

  private :: dp, i4

  integer(i4) :: Lx,Lt ! Lattice length size
  integer(i4) :: N_measurements
  integer(i4) :: N_thermalization
  integer(i4) :: N_skip
  character(20) :: algorithm
  logical :: equilibrium

  namelist /input_parameters/ Lx,Lt,N_thermalization, N_measurements, N_skip, algorithm, equilibrium

contains

  subroutine read_input_parameters()
    use iso_fortran_env, only : stdout => output_unit, stdin => input_unit, stderr => error_unit
    character(100) :: parameters_file
    integer(i4) :: inunit


    write(stdout,'(a)') "Enter the parameters file, please."
    read(stdin,'(a)') parameters_file
    write(stdout,'(a)') "User typed: ", trim(parameters_file)

    open(newunit = inunit, file = trim(parameters_file), status = 'old')
    read(inunit, nml = input_parameters)
    if( Lx <= 0) error stop "Lattice spatial length L must be > 0."
    if( Lt <= 0) error stop "Lattice time length L must be > 0."
    if( N_thermalization <= 0) error stop "Thermalization sweeps must be > 0."
    if( N_measurements <= 0) error stop "Number of measurements must be > 0."

    select case(algorithm)
    case("metropolis")
       print*, "ok"
    case("glauber")
       print*, "ok"
    case("heatbath")
       print*, "ok"
    case default
       error stop "Invalid algorithm"
    end select
    
    write(stdout, nml = input_parameters)


  end subroutine read_input_parameters

end module parameters
