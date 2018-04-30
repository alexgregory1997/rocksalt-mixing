module ionic_module
implicit none
contains
  subroutine distance(x,y,z,i,j,r,a)
    implicit none
    ! Local Variables
    integer, parameter :: dp=selected_real_kind(15,300)
    integer            :: x_distance, y_distance, z_distance
    ! Global Variables
    integer       :: x, y, z, i, j
    real(kind=dp) :: r, a

    ! Calculate distance in z direction.
    if (ceiling(real(i,dp)/(real(x,dp)*real(y,dp))).eq.ceiling(real(j,dp)/(real(x,dp)*real(y,dp)))) then
      z_distance=0
    else
      z_distance=ceiling(real(i,dp)/(real(x,dp)*real(y,dp)))-ceiling(real(j,dp)/(real(x,dp)*real(y,dp)))
    end if

    ! Calculate distance in x direction.
    if (mod(i,x).eq.mod(j,x)) then
      x_distance=0
    else
      if (mod(i,x).eq.0) then
        x_distance=x-mod(j,x)
      else if (mod(j,x).eq.0) then
        x_distance=mod(i,x)-x
      else
        x_distance=mod(i,x)-mod(j,x)
      end if
    end if

    ! Calculate distance in y direction.
    if (ceiling(real(i,dp)/real(y,dp)).eq.ceiling(real(j,dp)/real(y,dp))) then
      y_distance=0
    else if (mod(ceiling(real(i,dp)/real(y,dp)),z).eq.0) then
      y_distance=y-mod(ceiling(real(j,dp)/real(y,dp)),z)
    else if (mod(ceiling(real(j,dp)/real(y,dp)),z).eq.0) then
      y_distance=mod(ceiling(real(i,dp)/real(y,dp)),z)-y
    else
      y_distance=mod(ceiling(real(i,dp)/real(y,dp)),z)-mod(ceiling(real(j,dp)/real(y,dp)),z)
    end if

    ! Define the distance between points i and j.
    r=sqrt((real(z_distance,dp))**2+(real(x_distance,dp))**2+(real(y_distance,dp))**2)*a

  end subroutine distance

  ! This calculates the average bond length between ions within the lattice depending on the
  ! ratio of Mg to Ca ions present.
  subroutine lat_param(total_pos,a,lattice)
    implicit none
    ! Local Variables
    integer, parameter :: dp=selected_real_kind(15,300)
    integer            :: i, num_Mg, num_Ca
    real(kind=dp)      :: a_Mg, a_Ca
    ! Global Variables
    integer                            :: total_pos
    integer, dimension(:), allocatable :: lattice
    real(kind=dp)                      :: a

    ! Experimental values for MgO and CaO
    a_Mg=2.12e-10_dp
    a_Ca=2.4e-10_dp
    ! Initialise
    num_Mg=0
    num_Ca=0

    ! Calculate total number of Mg/Ca in the lattice
    do i=1,total_pos
      if (lattice(i).eq.-1) then
        num_Mg=num_Mg+1
      else if (lattice(i).eq.1) then
        num_Ca=num_Ca+1
      end if
    end do

    ! Calculate lattice paramter a for a given crystal
    a=((num_Mg*a_Mg)+(num_Ca*a_ca))/(num_Mg+num_Ca)

  end subroutine lat_param

  ! This calculates the interaction energy at a single point in the crystal
  subroutine charge(i,q,lattice)
    implicit none
    ! Local variables
    integer, parameter                 :: dp=selected_real_kind(15,300)
    real(kind=dp)                      :: q
    ! Global Variables
    real(kind=dp), parameter           :: e=1.60217662e-19_dp
    integer                            :: i
    integer, dimension(:), allocatable :: lattice

    ! Calculate charge on specific position in crystal
    ! Note - looks convoluted, but this is for possible extensions, i.e.
    ! looking at effect charge
    if (lattice(i).eq.0) then
      q=-2.0_dp*e
    else if (lattice(i).eq.-1) then
      q=2.0_dp*e
    else if (lattice(i).eq.1) then
      q=2.0_dp*e
    end if

  end subroutine charge

  ! This subroutine generates the repulsion parameters for a give pair of atoms.
  ! Rules are based off papers found - A uses geometric mean and row arithmetic.
  subroutine repulsion_param(A_MgMg,A_CaCa,A_OO,row_MgMg,row_CaCa,row_OO,A_param,row_param,pos,i,lattice)
    implicit none
    ! Local Variables
    integer, parameter :: dp=selected_real_kind(15,300)
    real(kind=dp)      :: A_pos, row_pos, A_on, row_on
    ! Global Variables
    real(kind=dp)                      :: A_MgMg, A_CaCa, A_OO, row_MgMg, row_CaCa, row_OO
    real(kind=dp)                      :: A_param, row_param
    integer                            :: pos, i
    integer, dimension(:), allocatable :: lattice

    if (lattice(pos).eq.-1) then
      A_pos=A_MgMg
      row_pos=row_MgMg
    else if (lattice(pos).eq.1) then
      A_pos=A_CaCa
      row_pos=row_CaCa
    else
      A_pos=A_OO
      row_pos=row_OO
    end if

    if (lattice(i).eq.-1) then
      A_on=A_MgMg
      row_on=row_MgMg
    else if (lattice(i).eq.1) then
      A_on=A_CaCa
      row_on=row_CaCa
    else
      A_on=A_OO
      row_on=row_OO
    end if

    A_param=sqrt(A_pos*A_on)
    row_param=0.5_dp*(row_pos+row_on)

  end subroutine repulsion_param

  ! This calculates the iterative part of the interaction energy
  function iterative_int(q_on,q_pos,r,A_param,row_param)
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    real(kind=dp), parameter :: pi=3.14159265359_dp, epsilon=8.85418782e-12_dp
    real(kind=dp)            :: q_on, q_pos
    real(kind=dp)      :: iterative_int, r, A_param, row_param

    iterative_int=(q_on/(4.0_dp*pi*epsilon))*(q_pos/r)+A_param*exp(-r/row_param)

  end function iterative_int

end module ionic_module
