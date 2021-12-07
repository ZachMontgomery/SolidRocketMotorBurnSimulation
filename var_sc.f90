module var_sc
implicit none
!***************************************************************************
integer, parameter :: sp = selected_real_kind(p=6)
integer, parameter :: wp = selected_real_kind(p=8)
!***************************************************************************
real(wp), parameter :: pi = 3.141592653589793_wp
! geometry variables
real(wp), parameter :: astar = 1.887_wp, exp_ratio = 3.746_wp, theta_exit =&
20._wp / 180._wp * pi, lport = 35._wp, rexit = sqrt(astar*exp_ratio/pi)
real(wp), parameter :: dout = 7.6_wp, rho_propellant = 1260._wp, din = &
3._wp, rthroat = sqrt(astar / pi), rg = 8314.4598_wp/23._wp
! combustion gas and burn properties
real(wp), parameter :: g = 1.18_wp, mol_weight = 23._wp, t0 = 2900._wp, n =&
0.16_wp, a = .12_wp, mach_crit = 0.15_wp, k = 1._wp, pamb = 86._wp
! rocket variables
real(wp), parameter :: lambda = .5_wp*(1._wp+cos(theta_exit))
!***************************************************************************
! variable port geometries
real(wp) :: rpo, rc, rfo, yo, acase, vcase
real(wp), dimension(0:4) :: s, rf, rp, y, theta, Ac, vport, aburn, afuel, &
vfuel, mass_propellant, breach_time, burn_time
integer :: cracks
!***************************************************************************
! integration variables
integer :: cont
real(wp) :: dt, time
real(wp), dimension(0:4,2) :: x
! Index 0:4 refers to the scenario being run. The index is the number of
! cracks. The other column is the state variables Po and r
!***************************************************************************
! plot variables
integer, parameter :: res = 1000
real, dimension(res,0:4,2) :: outline
real, dimension(res) :: casex, casey
! res reperesent number of point to make the outline
! 0:4 reperensent which scenario the array is for, how many cracks
! 2 reperesents 1 for r and 2 for theta
real, allocatable, dimension(:,:) :: chamber_pressure, regression_rate, &
massflow, thrust, isp, area2volume, isp_mean
real, allocatable, dimension(:) :: timep
integer :: siz
contains
!***************************************************************************
function exp_ratio2mach(g, rt, r, guess, error, max_iterations)
! returns mach given radius at the throat, gamma and radius at desired point
integer, parameter :: prec = selected_real_kind(8)
real(prec), intent(in) :: g, rt, r, guess, error
integer, intent(in) :: max_iterations
real(prec) :: f, df, exp_ratio2mach, a, b, c, d, xo, xn, ea
integer :: counter
ea = 1._prec
counter = 0
xo = guess
a = r**2 / rt**2
b = (g+1._prec)/2._prec/(g-1._prec)
c = (1._prec - 3._prec * g) / (2._prec - 2._prec * g)
do while(ea > error)
	if(counter >= max_iterations) then
		write(*,*)
		write(*,*) 'Newtons Method inside of the exp_ratio2mach Function'
		write(*,*) 'Failed to Converge within ',counter,' iterations.'
		write(*,*) 'Approximate error was ',ea
		write(*,*)
		read(*,*)
	end if
	d = xo**2
	f = (2._prec/(g+1._prec)*(1._prec+(g-1._prec)/2._prec*d))**b/xo-a
	df = 2._prec**c*(d-1._prec)/d/(2._prec+d*(g-1._prec))*((1._prec+(g-&
			1._prec)/2._prec*d)/(g+1._prec))**b
	xn = xo - f / df
	ea = abs(xn - xo) / xo
	xo = xn
	counter = counter + 1
end do
exp_ratio2mach = xn
end function exp_ratio2mach
!***************************************************************************
end module var_sc
