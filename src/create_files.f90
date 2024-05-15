module create_files
  use check_files_directories_mod
  use number2string_mod
  implicit none

contains

  subroutine create_measurements_file(Lx1,Lt1,beta,algorithm1,equilibrium1)
    use parameters
    integer, intent(in) :: Lx1,Lt1
    real(8), intent(in) :: beta
    character(*), intent(in) :: algorithm1
    logical :: equilibrium1, file_exists, condition
    integer :: outunit, i

    character(100) :: directory, eq, data_file

    if(equilibrium1 .eqv. .true.) then
       eq = "equilibrium"
    else
       eq = "out_of_equilibrium"
    end if
    
    directory = "data"
    call check_directory(trim(directory))
    
    directory = trim(directory)//"/Lx="//trim(int2str(Lx1))
    call check_directory(trim(directory))
    
    directory = trim(directory)//"/Lt="//trim(int2str(Lt1))
    call check_directory(trim(directory))
    
    directory = trim(directory)//"/"//trim(eq)
    call check_directory(trim(directory))

    directory = trim(directory)//"/"//trim(algorithm1)
    call check_directory(trim(directory))

    directory = trim(directory)//"/beta="//real2str(beta)
    call check_directory(trim(directory))

    directory = trim(directory)//"/observables_"

    i = 1
    do 
       data_file = trim(directory)//trim(int2str(i))//".dat"
       call check_file(trim(data_file), file_exists)
       condition = file_exists
       if (condition .eqv. .false.) exit
       i = i + 1
    end do
    open(unit = 100, file = data_file)
    write(100, nml = input_parameters )
    
  end subroutine create_measurements_file

end module create_files
