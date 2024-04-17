module check_files_directories_mod

    implicit none

    contains

    subroutine check_directory(directory)

    implicit none
        character(*), intent(in) :: directory
        logical ::  dir_exists


        inquire(file=directory, exist=dir_exists) ! ask wether the directory exists or not

        if(dir_exists .eqv. .false.) call execute_command_line('mkdir '//directory) ! If not it creates it


    end subroutine check_directory

    subroutine check_file(filepath,file_exists)

        character(*), intent(in) :: filepath
        logical, intent(out) :: file_exists

        inquire(file = filepath, exist = file_exists)
        if(file_exists .eqv. .false.)then
           call execute_command_line('touch '//filepath)
           !file_exists = .true.
        end if
    end subroutine check_file


end module check_files_directories_mod
