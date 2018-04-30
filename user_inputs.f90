module user_inputs
implicit none
contains
  subroutine initialise_lattice(total_pos,form,lattice)
    implicit none
    integer :: i, total_pos
    character(len=6) :: form
    integer, dimension(:), allocatable :: lattice

    print*, 'Do you wish for the initial plane to contain five cations or anions?'
    print*, 'Please enter cation or anion: '
    read*, form
    select case (form)
    case ('cation')
      do i=1,total_pos
        if (mod(i,2).eq.0) then
          lattice(i)=0
        else
          lattice(i)=-1  ! Need to add something to select Ca/Mg
        end if
      end do
    case ('anion')
      do i=1,total_pos
        if (mod(i,2).eq.0) then
          lattice(i)=1
        else
          lattice(i)=0
        end if
      end do
    case default
      print*, 'Program terminated. Invalid entry. Please try again.'
    end select
  end subroutine initialise_lattice

end module user_inputs
