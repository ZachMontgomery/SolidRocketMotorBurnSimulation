program solid_crack
use var_sc
use dislin
implicit none
!***************************************************************************
! declare variables
character :: live_plot
!***************************************************************************
! initialize variables
x(:,1) = pamb
x(:,2) = .5_wp * din
rc = .5_wp * dout
rpo = .5_wp * din
yo = .04_wp
rfo = 3._wp
acase = pi * rc**2
vcase = pi * rc**2 * lport
do cracks = 0, 4
	call update_geometry(cracks)
end do
burn_time = 0._wp
breach_time = 0._wp
time = 0._wp
! get variables from user
write(*,*)
write(*,*) 'Enter value for time step: '
write(*,*)
read(*,*) dt
write(*,*)
write(*,*) 'Watch live plot (y/n): '
write(*,*)
read(*,*) live_plot
! allocate variables and initialize
siz = nint(10._wp/dt)
allocate(chamber_pressure(siz,0:4), regression_rate(siz,0:4), massflow(siz,&
0:4), thrust(siz,0:4), isp(siz,0:4))
allocate(area2volume(siz,0:4), isp_mean(siz,0:4), timep(siz))
chamber_pressure = pamb
regression_rate = 0.
massflow = 0.
thrust = 0.
isp = 0.
do cracks = 0, 4
	area2volume(:,cracks) = aburn(cracks) / vport(cracks)
end do
isp_mean = 0.
timep = 0.
!***************************************************************************
call metafl('xwin')
!call scrmod('reverse')
call window(25,25,2000,1000)
call page(7000,3500)
call disini
call bufmod('off','sendbf')
do cracks = 0, res-1
	casex(cracks+1)=real(rc,sp)*cos(2.*real(pi,sp)/real(res-1)*real(cracks))
	casey(cracks+1)=real(rc,sp)*sin(2.*real(pi,sp)/real(res-1)*real(cracks))
end do
! integrate and plot
cont = 0
siz = 0
do while (cont < 5)
	siz = siz + 1
	cont = 0
	!write(*,*) 'Time = ',time
	do cracks = 0, 4
		if (mass_propellant(cracks) <= 0._wp) then
			cont = cont + 1
			if (burn_time(cracks) == 0._wp) burn_time(cracks) = time
		else
			call trapezoidal(cracks)
			call update_geometry(cracks)
			call perimeter(cracks)
			call rocket_calcs(siz,cracks)
		end if
	end do
	time = time + dt
	timep(siz) = time
	! plots
	if (live_plot == 'y') then
		call plotter
	else
		if (mod(real(time,sp),1.) == 0.) then
			call erase
			call number(real(time,sp),3,3500,1750)
			call sendbf
		end if
	end if
end do
call plotter
!***************************************************************************
call disfin
!***************************************************************************
write(*,*) 'Breach times are'
write(*,*) breach_time
write(*,*) 'Total burn times are'
write(*,*) burn_time
!read(*,*)
!***************************************************************************
deallocate(chamber_pressure, regression_rate, massflow, thrust, isp, &
area2volume, isp_mean, timep)
!***************************************************************************
end program solid_crack
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
!***************************************************************************
subroutine update_geometry(i)
use var_sc
integer, intent(in) :: i
s(i) = x(i,2) - rpo
if (rf(i) < rc) rf(i) = rfo + s(i)
if (rp(i) < rc) rp(i) = rpo + s(i)
y(i) = yo + s(i)
theta(i) = atan(2._wp * y(i) / (rf(i) + rp(i)))
ac(i) = pi * rp(i)**2 + real(i)*theta(i)*(rf(i)**2 - rp(i)**2)
vport(i) = ac(i) * lport
if (rf(i) < rc) then
	aburn(i) = lport * 2._wp * (pi * rp(i) + (rf(i) - rp(i) + theta(i) * &
				(rf(i) - rp(i))) * real(i))
else
	aburn(i) = lport * 2._wp * (pi * rp(i) + (rf(i) - rp(i) - theta(i) * &
				rp(i)) * real(i))
	if (breach_time(i) == 0._wp) breach_time(i) = time
end if
afuel(i) = acase - ac(i)
vfuel(i) = afuel(i)*lport
mass_propellant(i) = rho_propellant * vfuel(i) / 100._wp**3
!write(*,*) 'mass of ',i,' is ',mass_propellant(i)
!read(*,*)
end subroutine update_geometry
!***************************************************************************
subroutine trapezoidal(i)
use var_sc
use func
integer, intent(in) :: i
real(wp), dimension(2) :: xhat
!write(*,*) 'inside trapezoidal'
!write(*,*) 'x = ',x(i,:)
!read(*,*)
!write(*,*) 'before derivative'
xhat = x(i,:) + dt * derivative(x(i,:),i)
!write(*,*) 'after derivative'
!write(*,*) 'xhat = ',xhat
x(i,:) = x(i,:)+dt/2._wp*(derivative(x(i,:),i)+derivative(xhat,i))
!write(*,*) 'new x = ',x(i,:)
!read(*,*)
end subroutine trapezoidal
!***************************************************************************
subroutine perimeter(i)
use var_sc
integer, intent(in) :: i
real(wp) :: t
integer :: j
do j = 0, res-1
	! calculate theta
	t = 2._wp * pi / real(res-1,wp) * real(j,wp)
	outline(j+1,i,2) = t
	! calculate r
	if (i > 0 .and. (t <= theta(i) .or. t >= 2._wp*pi - theta(i))) then
		outline(j+1,i,1) = rf(i)
	else if (i > 1 .and. (t <= pi + theta(i) .and. t >= pi - theta(i))) then
		outline(j+1,i,1) = rf(i)
	else if (i>2 .and. t<=pi/2._wp+theta(i) .and. t>=pi/2._wp-theta(i)) then
		outline(j+1,i,1) = rf(i)
	else if (i>3 .and. t<=3._wp*pi/2._wp+theta(i) .and. t>=3._wp*pi/2._wp-&
															theta(i)) then
		outline(j+1,i,1) = rf(i)
	else
		outline(j+1,i,1) = rp(i)
	end if
end do
end subroutine perimeter
!***************************************************************************
subroutine plotter
use var_sc
use func
use dislin
real :: u, v
integer :: e
u = 1.1*rc
v = rc / 5.
call erase
! no cracks
e = 0
call color('fore')
call axslen(1000,1000)
call axspos(250,1250)
call labels('float','y')
call labtyp('vert','x')
call labdig(1,'y')
call name('','xy')
call graf(-u,u,-u,v,-u,u,-u,v)
call curve(casex,casey,res)
call color('red')
call curve(pol2cart(outline(:,e,1),outline(:,e,2),res,'x'),&
			pol2cart(outline(:,e,1),outline(:,e,2),res,'y'),res)
call endgrf
! 1 crack
e = 1
call color('fore')
call axslen(1000,1000)
call axspos(1250,1250)
call labels('none','y')
call labtyp('vert','x')
call graf(-u,u,-u,v,-u,u,-u,v)
call curve(casex,casey,res)
call color('orange')
call curve(pol2cart(outline(:,e,1),outline(:,e,2),res,'x'),&
			pol2cart(outline(:,e,1),outline(:,e,2),res,'y'),res)
call endgrf
! 2 cracks
e = 2
call color('fore')
call axslen(1000,1000)
call axspos(2250,1250)
call labtyp('vert','x')
call graf(-u,u,-u,v,-u,u,-u,v)
call curve(casex,casey,res)
call color('yellow')
call curve(pol2cart(outline(:,e,1),outline(:,e,2),res,'x'),&
			pol2cart(outline(:,e,1),outline(:,e,2),res,'y'),res)
call endgrf
! 3 cracks
e = 3
call color('fore')
call axslen(1000,1000)
call axspos(3250,1250)
call labtyp('vert','x')
call graf(-u,u,-u,v,-u,u,-u,v)
call curve(casex,casey,res)
call color('green')
call curve(pol2cart(outline(:,e,1),outline(:,e,2),res,'x'),&
			pol2cart(outline(:,e,1),outline(:,e,2),res,'y'),res)
call endgrf
! 4 cracks
e = 4
call color('fore')
call axslen(1000,1000)
call axspos(4250,1250)
call labtyp('vert','x')
call graf(-u,u,-u,v,-u,u,-u,v)
call curve(casex,casey,res)
call color('blue')
call curve(pol2cart(outline(:,e,1),outline(:,e,2),res,'x'),&
			pol2cart(outline(:,e,1),outline(:,e,2),res,'y'),res)
call endgrf
! chamber pressure
call color('fore')
call axslen(1000,1000)
call axspos(400,2500)
call labels('float','y')
call labtyp('vert','x')
call name('Time (s)','x')
call name('Chamber Pressure (kPa)','y')
v = 0.
do e = 0, 4
	u = maxval(chamber_pressure(:,e))
	if (u > v) v = u
end do
u = time
call graf(0.,u,0.,.5,0.,1.1*v,0.,.1*v)
call color('red')
call curve(timep(1:siz),chamber_pressure(1:siz,0),siz)
call color('orange')
call curve(timep(1:siz),chamber_pressure(1:siz,1),siz)
call color('yellow')
call curve(timep(1:siz),chamber_pressure(1:siz,2),siz)
call color('green')
call curve(timep(1:siz),chamber_pressure(1:siz,3),siz)
call color('blue')
call curve(timep(1:siz),chamber_pressure(1:siz,4),siz)
call endgrf
! regression rate
call color('fore')
call axslen(1000,1000)
call axspos(1700,2500)
call labtyp('vert','x')
call labdig(2,'y')
call name('Regression Rate (cm/s)','y')
v = 0.
do e = 0, 4
	u = maxval(regression_rate(:,e))
	if (u > v) v = u
end do
u = time
call graf(0.,u,0.,.5,0.,1.1*v,0.,.1*v)
call color('red')
call curve(timep(1:siz),regression_rate(1:siz,0),siz)
call color('orange')
call curve(timep(1:siz),regression_rate(1:siz,1),siz)
call color('yellow')
call curve(timep(1:siz),regression_rate(1:siz,2),siz)
call color('green')
call curve(timep(1:siz),regression_rate(1:siz,3),siz)
call color('blue')
call curve(timep(1:siz),regression_rate(1:siz,4),siz)
call endgrf
! massflow
call color('fore')
call axslen(1000,1000)
call axspos(3000,2500)
call labtyp('vert','x')
call name('Massflow (kg/s)','y')
v = 0.
do e = 0, 4
	u = maxval(massflow(:,e))
	if (u > v) v = u
end do
u = time
call graf(0.,u,0.,.5,0.,1.1*v,0.,.1*v)
call color('red')
call curve(timep(1:siz),massflow(1:siz,0),siz)
call color('orange')
call curve(timep(1:siz),massflow(1:siz,1),siz)
call color('yellow')
call curve(timep(1:siz),massflow(1:siz,2),siz)
call color('green')
call curve(timep(1:siz),massflow(1:siz,3),siz)
call color('blue')
call curve(timep(1:siz),massflow(1:siz,4),siz)
call endgrf
! thrust
call color('fore')
call axslen(1000,1000)
call axspos(4400,2500)
call labtyp('vert','x')
call labdig(1,'y')
call name('Thrust (N)','y')
v = 0.
do e = 0, 4
	u = maxval(thrust(:,e))
	if (u > v) v = u
end do
u = time
call graf(0.,u,0.,.5,0.,1.1*v,0.,.1*v)
call color('red')
call curve(timep(1:siz),thrust(1:siz,0),siz)
call color('orange')
call curve(timep(1:siz),thrust(1:siz,1),siz)
call color('yellow')
call curve(timep(1:siz),thrust(1:siz,2),siz)
call color('green')
call curve(timep(1:siz),thrust(1:siz,3),siz)
call color('blue')
call curve(timep(1:siz),thrust(1:siz,4),siz)
call endgrf
! send plot data to screen
call sendbf
end subroutine plotter
!***************************************************************************
subroutine rocket_calcs(i,j)
use var_sc
use func
integer, intent(in) :: i,j
real, dimension(2) :: dumb
real(wp) :: pexit, mexit, texit
! calculate chamber pressure
chamber_pressure(i,j) = x(j,1)
! calculate regression rate
dumb = derivative(x(j,:),j)
regression_rate(i,j) = dumb(2)
! calculate massflow
massflow(i,j) = rho_propellant * aburn(j) * dumb(2) / 100.**3
! calculate thrust
mexit = exp_ratio2mach(g,rthroat,rexit,1.1_wp,.000001_wp,500)
texit = t0/(1._wp+(g-1._wp)/2._wp*mexit**2)
pexit = x(j,1)/(1._wp+(g-1._wp)/2._wp*mexit**2)**(g/(g-1._wp))
thrust(i,j) = lambda*massflow(i,j)*mexit*sqrt(g*rg*texit)+astar*exp_ratio*(&
			pexit-pamb)/1000._wp
! calculate isp
!isp(i,j) = 
end subroutine rocket_calcs
!***************************************************************************
