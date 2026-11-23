! Test-only wrapper for ISORROPIA II in CMAQ 5.5.0.4, aero6/isocom.f.
subroutine mkpp_reference_cmaq_ii(totals, rh, temperature, gas, particle, sulfate, water, hydronium, status) bind(C)
  use iso_c_binding, only: c_double, c_int
  implicit none
  real(c_double), intent(in) :: totals(8)
  real(c_double), value :: rh, temperature
  real(c_double), intent(out) :: gas(3), particle(2), sulfate, water, hydronium
  integer(c_int), intent(out) :: status
  real(c_double) :: control(2), total_out(8), liquid(15), solid(19), diagnostics(9)
  character(len=15) :: regime
  character(len=40) :: messages(25)
  integer :: errors(25), nerrors
  logical :: overflow

  control = [0.0_c_double, 1.0_c_double] ! Forward, metastable
  ! Use tabulated activities, same numerical controls as original Lite, mass adjustment off.
  call setparm(2, 1, 1.0d-6, 100, 4, 5.0d-2, 5, 0)
  call isoropia(totals, rh, temperature, control, total_out, gas, liquid, solid, regime, diagnostics)
  call iserrinf(errors, messages, nerrors, overflow)
  status = nerrors
  if (overflow) status = -1

  ! These cases contain only NHx, sulfate and nitrate. Count all salts of those components.
  particle(1) = liquid(3) + liquid(9) + solid(2) + solid(4) + 2*solid(6) + solid(8) + 3*solid(9)
  particle(2) = liquid(7) + liquid(11) + solid(1) + solid(2) + 2*solid(11) + solid(15) + 2*solid(18)
  sulfate = liquid(5) + liquid(6) + sum(solid(5:8)) + 2*solid(9) + solid(10) + solid(13) + solid(14) + solid(17)
  water = 0.018_c_double * liquid(8) ! mol H2O/m3 to kg/m3
  hydronium = liquid(1)              ! mol H+/m3 air
end subroutine mkpp_reference_cmaq_ii
