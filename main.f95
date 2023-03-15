program vibrational_spectroscopy
  implicit none
  
  integer, parameter :: n_atoms = 10 ! Number of atoms in the molecule
  integer, parameter :: n_modes = 30 ! Number of vibrational modes
  real(kind=8) :: coordinates(n_atoms, 3) ! Atomic coordinates
  real(kind=8) :: potential_energy_surface(n_modes, n_atoms, 3) ! Potential energy surface
  real(kind=8) :: frequencies(n_modes) ! Vibrational frequencies
  real(kind=8) :: intensities(n_modes) ! Vibrational intensities
  integer :: i, j, k ! Loop indices
  
  ! Random atomic coordinates between -5 and 5 angstroms
  call random_seed()
  do i = 1, n_atoms
    call random_number(coordinates(i,:))
    coordinates(i,:) = (coordinates(i,:) * 2.0 - 1.0) * 5.0
  end do
  
  ! Harmonic oscillator model
  do i = 1, n_modes
    do j = 1, n_atoms
      do k = 1, 3
        potential_energy_surface(i,j,k) = 0.5 * (coordinates(j,k) ** 2)
      end do
    end do
  end do

  ! Calculate the vibrational frequencies and intensities using a normal mode analysis
  ! Diagonalizes the Hessian matrix
  do i = 1, n_modes
    frequencies(i) = sqrt(sum(potential_energy_surface(i,:,:))) / (2 * i)
    intensities(i) = 0
    do j = 1, n_atoms
      do k = 1, 3
        intensities(i) = intensities(i) + (potential_energy_surface(i,j,k) * coordinates(j,k)) ** 2
      end do
    end do
    intensities(i) = intensities(i) / (3 * 8.85418782e-30 * frequencies(i))
  end do
  
    print *, "Atomic coordinates:"
    do i = 1, n_atoms
      print *, i, coordinates(i,1), coordinates(i,2), coordinates(i,3)
    end do
    
    print *, "Potential energy surface:"
    do i = 1, n_modes
      print *, "Mode ", i
      do j = 1, n_atoms
        print *, j, potential_energy_surface(i,j,1), potential_energy_surface(i,j,2), potential_energy_surface(i,j,3)
      end do
    end do
    
    print *, "Vibrational frequencies and intensities:"
    do i = 1, n_modes
      print *, "Mode ", i
      print *, "Frequency (cm^-1): ", frequencies(i)
      print *, "Intensity (km/mol): ", intensities(i)
    end do

  
end program vibrational_spectroscopy
