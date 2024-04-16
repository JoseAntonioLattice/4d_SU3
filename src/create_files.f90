module create_files
  use check_files_directories_mod
  use number2string_mod
  implicit none

contains

  subroutine create_measurements_file(Lx,Lt,beta,algorithm,equilibrium)
    integer, intent(in) :: Lx,Lt
    real(8), intent(in) :: beta
    character(*), intent(in) :: algorithm
    logical :: equilibrium
    integer :: outunit

    character(100) :: directory, eq

    if(equilibrium .eqv. .true.) then
       eq = "equilibrium"
    else
       eq = "out_of_equilibrium"
    end if
    
    directory = "data"
    call check_directory(trim(directory))
    
    directory = trim(directory)//"/Lx="//trim(int2str(Lx))
    call check_directory(trim(directory))

    
    directory = trim(directory)//"/Lt="//trim(int2str(Lt))
    call check_directory(trim(directory))
    
    directory = trim(directory)//"/"//trim(eq)
    call check_directory(trim(directory))

    directory = trim(directory)//"/"//trim(algorithm)
    call check_directory(trim(directory))

    directory = trim(directory)//"/beta="//real2str(beta)
    call check_directory(trim(directory))

    
    open(unit = 100, file = trim(directory)//'/observables.dat')
    
  end subroutine create_measurements_file

end module create_files
