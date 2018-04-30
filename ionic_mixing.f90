program ionic_mixing
  use ionic_module
  use user_inputs
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)
  real(kind=dp), parameter :: pi=3.14159265359_dp, epsilon=8.85418782e-12_dp, e=1.60217662e-19_dp, k=1.38064852e-23_dp
  integer            :: ierr, ios
  integer            :: total_pos, x, y, z, i, j, m, ms, max
  integer            :: on
  real(kind=dp)      :: q_pos, q_on
  real(kind=dp)      :: r, a, initial_energy, final_energy, pair_energy, total_energy, beta, T, ran1, ran2
  real(kind=dp)      :: A_param, row_param
  real(kind=dp)      :: A_MgMg, A_CaCa, A_OO, row_MgMg, row_CaCa, row_OO
  integer, dimension(:), allocatable :: lattice
  character(len=6)   :: form

  ! Initialise Variables
  ! Number of lattice positions in length, width and height.
  x=15
  y=15
  z=15
  ! Total number of lattice postions within crystal.
  total_pos=x*y*z
  ! Distance between two points
  r=0.0_dp  !probably don't need to initialise
  ! Repulsion coefficents - find experimental values
  A_OO=2.76e-10_dp ! from paper
  A_CaCa=2.0e-10_dp ! from paper
  A_MgMg=1.44e-10_dp! from paper
  row_OO=0.2192e-10_dp ! from paper
  row_CaCa=0.08e-10_dp ! from paper
  row_MgMg=0.1e-10_dp ! yet to find
  ! Temperature dependent factor
  T=1.0_dp  ! Kelvin
  beta=1/(k*T)
  ! Monte Carlo Step - when to calculate total lattice energy
  ms=50
  ! How many iterations
  max=1000

  ! Allocate Arrays
  ! Allocate lattice as a one-dimensional chain with size 'total_pos', defined above.
  allocate(lattice(total_pos),stat=ierr)
  if (ierr.ne.0) stop 'error allocating lattice array'

  ! Initialise Arrays
  call initialise_lattice(total_pos,form,lattice)

  ! Calculates lattice paramter for a given system
  call lat_param(total_pos,a,lattice)

  open(unit=50,file='test.dat',iostat=ios)
  if (ios.ne.0) stop 'error opening test.dat'

  ! Monte Carlo Process
  do m=1,max

    ! Select random atom that is not oxygen
    do
      call random_number(ran1)
      on=int(total_pos*ran1)
      if (on.lt.1 .or. on.gt.total_pos) then
        cycle
      else if (lattice(on).eq.0) then
        cycle
      else
        exit
      end if
    end do

    ! Calculate charge at this point and lattice parameter for this system
    call charge(on,q_on,lattice)
    call lat_param(total_pos,a,lattice)

    ! Calculate interaction energy at this point
    initial_energy=0.0_dp
    do i=1,total_pos
      if (i.eq.on) then
        initial_energy=initial_energy+0.0_dp
      else
        call charge(i,q_pos,lattice)
        call distance(x,y,z,on,i,r,a)
        call repulsion_param(A_MgMg,A_CaCa,A_OO,row_MgMg,row_CaCa,row_OO,A_param,row_param,on,i,lattice)
        initial_energy=initial_energy+iterative_int(q_on,q_pos,r,A_param,row_param)
        print*, initial_energy, r
      end if
    end do

    ! Change atom type at this point
    lattice(on)=lattice(on)*(-1)

    ! Recalculate charge at this point and lattice paramter for the system
    call charge(on,q_on,lattice)
    call lat_param(total_pos,a,lattice)

    ! Recalculate interaction energy at this point
    final_energy=0.0_dp
    do i=1,total_pos
      if (i.eq.on) then
        final_energy=final_energy+0.0_dp
      else
        call charge(i,q_pos,lattice)
        call distance(x,y,z,on,i,r,a)
        call repulsion_param(A_MgMg,A_CaCa,A_OO,row_MgMg,row_CaCa,row_OO,A_param,row_param,on,i,lattice)
        final_energy=final_energy+iterative_int(q_on,q_pos,r,A_param,row_param)
      end if
    end do

    ! Monte Carlo algorithm - decide whether to keep change or not
    if (final_energy-initial_energy.gt.0) then
      call random_number(ran2)
      if (exp(-(final_energy-initial_energy)*beta).lt.ran2) then
        print*, exp(-(final_energy-initial_energy)*beta), ran2
        lattice(on)=lattice(on)*(-1)
    end if

    !print*, on, q_on, a, initial_energy, final_energy, final_energy-initial_energy, row_param

    ! Calculating total energy of the system every ms amount of changes
    if (mod(m,ms).eq.0) then
      call lat_param(total_pos,a,lattice)
      total_energy=0.0_dp
      do i=1,total_pos
        call charge(i,q_on,lattice)
        pair_energy=0.0_dp
        do j=1,total_pos
          if (j.eq.i) then
            pair_energy=pair_energy+0.0_dp
          else
            call charge(j,q_pos,lattice)
            call distance(x,y,z,i,j,r,a)
            call repulsion_param(A_MgMg,A_CaCa,A_OO,row_MgMg,row_CaCa,row_OO,A_param,row_param,i,j,lattice)
            pair_energy=pair_energy+iterative_int(q_on,q_pos,r,A_param,row_param)
          end if
        end do
      total_energy=total_energy+(pair_energy/2.0_dp)
      end do
      write(50,*) m, total_energy
    end if
  end do

  close(50,iostat=ios)
  if (ios.ne.0) stop 'error closing test.dat'

  if (sum(lattice).eq.0) then
    print*, 'Equal number of Ca and Mg'
  else if (sum(lattice).gt.0) then
    print*, 'More Ca atoms to Mg'
    print*, 'Lattice sum: ', sum(lattice)
  else
    print*, 'More Mg atoms to Ca'
    print*, 'Lattice sum: ', sum(lattice)
  end if

  ! Figure a way out to plot this as 3d - may need a function to work out coordinates
  open(unit=50,file='3d_data.dat',iostat=ios)
  if (ios.ne.0) stop 'error opening 3d_data.dat'

  do i=1,total_pos
    write(50,*) lattice(i)
  end do

  close(50,iostat=ios)
  if (ios.ne.0) stop 'error closing 3d_data.dat'

  print*, 'FINISHED'

end program ionic_mixing
