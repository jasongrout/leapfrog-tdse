program leapfrog
  implicit none

  integer :: nx, nt, i, t
  real :: xmin, xmax, x, dx
  real :: tmin, tmax, dt
  real :: k0, w, xc, s, norm

  complex :: wave
  real, allocatable :: psi_r(:,:), psi_im(:,:)
  complex, external :: wavefn

  ! wavepacket parameters
  k0 = 0.0
  w = 0.05
  xc = 1.5

  ! grid parameters
  xmin = 1.0
  xmax = 2.0
  tmin = 0.0
  tmax = 1e-2
  nt = 1e4
  nx = 500
  dx = (xmax - xmin) / nx
  dt = (tmax - tmin) / nt

  s = dt / dx ** 2

  if (s .gt. 0.5) then
     ! 0.5 was just an arbitrary value mentioned in some of the literature
     ! Actually calculating the correct convergence threshold is difficult.
     write(*,*), "dt / dx^2 is off, unstable solution"
     stop
  end if

  allocate(psi_r(nx, nt), psi_im(nx, nt))

  ! Store initial wavefunction
  do i = 1, nx
     x = xmin + i * dx
     wave = wavefn(x, xc, w, k0)
     psi_r(i, 1) = real(wave)
     psi_r(i, 2) = psi_r(i, 1)
     psi_im(i, 1) = aimag(wave)
     psi_im(i, 2) = psi_im(i, 1)
  end do

  ! Leapfrog method, described at http://helium.bradley.edu/PICUP/doc/tdse/Numerical_Strategy.pdf
  ! To do - right now, psi_r and psi_im are twice as large as they need to be
  do t = 3, nt
     do i = 1, nx
        if (mod(t, 2) .eq. 1) then
           ! Update real part at even time steps
           psi_r(i, t) = psi_r(i, t - 2) - s * psi_im(i + 1, t - 1) - s * psi_im(i - 1, t - 1) + 2 * s * psi_im(i, t - 1)
        else
           ! Update complex part at odd time steps
           psi_im(i, t) = psi_im(i, t - 2) + s * psi_r(i + 1, t - 1) + s * psi_r(i - 1, t - 1) - 2 * s * psi_r(i, t - 1)
        end if
     end do
  end do

  open(3, file = "leapfrog_norm.out")

  ! write output
  ! to do - right now, we just write out the norm at each point at each double time step. This is sub-optimal.
  do i = 1, nx
     x = xmin + i * dx
     write(3, "(1f8.4)", advance = "no"), x
     do t = 1, nt / 2-1
        ! We can't use the typical norm here
        if (mod(t, 2) .eq. 0) then
           norm = psi_r(i, 2 * t - 1) ** 2 + psi_im(i, 2 * (t - 1)) * psi_im(i, 2 * (t + 1))
        else
           norm = psi_im(i, 2 * t) ** 2 + psi_r(i, 2 * (t - 1) - 1) * psi_r(i, 2 * (t + 1) - 1)
        end if
        write(3, "(1f8.4)", advance="no"), norm
     end do
     write(3, *)
  end do

  close(3)
  deallocate(psi_r, psi_im)

end program leapfrog

function wavefn(x, xc, w, k0) result(res)
  ! Returns the value of a Gaussian packet at a specific point.
  implicit none
  real, intent(in) :: x, xc, w, k0
  complex :: res

  res = cmplx(-(x - xc) ** 2 / w ** 2, k0 * x)
  res = exp(res)
end function wavefn
