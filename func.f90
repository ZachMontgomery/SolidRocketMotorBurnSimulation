module func
use var_sc

contains
!***************************************************************************
function derivative(state,i)
use var_sc
integer, intent(in) :: i
real(wp), intent(in), dimension(2) :: state
real(wp), dimension(2) :: derivative
real(wp) :: Po, Podot, rdot, mach_port, dummy
Po = state(1)
mach_port = exp_ratio2mach(g,rthroat,state(2),0.01_wp,.000001_wp,500)
if (mach_port < 0._wp) then
	write(*,*) 'mach port = ',mach_port
	read(*,*)
end if
rdot = a * Po ** n * (1._wp + k * mach_port / mach_crit) / (1._wp + k)
dummy = aburn(i)*rdot / vport(i) * (rho_propellant * rg * t0 / 1000._wp - Po)
Podot =  dummy - astar / vport(i) * Po * sqrt(g * rg * t0 * (2._wp / (g + &
		1._wp))**((g + 1._wp) / (g - 1._wp))) * 100._wp
derivative(1) = Podot
derivative(2) = rdot
end function
!***************************************************************************
function pol2cart(r,angle,length,xory)
integer, intent(in) :: length
character, intent(in) :: xory
real, dimension(length), intent(in) :: r, angle
real, dimension(length) :: pol2cart
integer :: counter
if (xory == 'x') then
	do counter = 1, length
		pol2cart(counter) = r(counter) * cos(angle(counter))
	end do
else if (xory == 'y') then
	do counter = 1, length
		pol2cart(counter) = r(counter) * sin(angle(counter))
	end do
else
	write(*,*) 'Wrong choice in pol2cart function'
	read(*,*)
end if
end function pol2cart
!***************************************************************************
end module func
