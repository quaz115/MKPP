! Test-only wrapper for the original ISORROPIA-Lite fixed-form Fortran.
! ISOROPIA's arrays and units are documented in isolite1_0_src/isocom.f.
subroutine mkpp_reference_lite(totals, organic, rh, temperature, gases, particles, sulfate, water, status) bind(C)
  use iso_c_binding, only: c_double, c_int
  implicit none
  real(c_double), intent(in) :: totals(8), organic(3)
  real(c_double), value :: rh, temperature
  real(c_double), intent(out) :: gases(3), particles(3), sulfate, water(2)
  integer(c_int), intent(out) :: status
  real(c_double) :: control(1), total_out(8), liquid(39), solid(19), diagnostics(9)
  character(len=15) :: regime
  character(len=40) :: messages(25)
  integer :: errors(25), nerrors
  logical :: overflow

  control = 0.0_c_double  ! forward, metastable aerosol
  ! Original distributed isrpia.cnf settings, including mass adjustment off.
  call setparm(2, 1.0d-6, 100, 4, 5.0d-2, 5, 0)
  call isoropia(totals, organic, rh, temperature, control, total_out, gases, liquid, solid, regime, diagnostics)
  call iserrinf(errors, messages, nerrors, overflow)
  status = nerrors
  if (overflow) status = -1

  ! Accumulate full condensed families independently of a total-minus-gas calculation.
  ! The first term is the ion, followed by undissociated dissolved gas and solids.
  particles(1) = liquid(3) + liquid(9) + solid(2) + solid(4) + 2*solid(6) + solid(8) + 3*solid(9)
  particles(2) = liquid(7) + liquid(11) + solid(1) + solid(2) + 2*solid(11) + solid(15) + 2*solid(18)
  particles(3) = liquid(4) + liquid(10) + solid(3) + solid(4) + 2*solid(12) + solid(16) + 2*solid(19)
  ! Sulfate and bisulfate ions, plus all sulfate-bearing salts. Letovicite has two sulfates.
  sulfate = liquid(5) + liquid(6) + sum(solid(5:8)) + 2*solid(9) + solid(10) + &
            solid(13) + solid(14) + solid(17)
  ! The legacy result expresses water in mol/m3. The C++ host reports kg/m3.
  water = 0.018_c_double * [liquid(8), liquid(39)]
end subroutine mkpp_reference_lite
