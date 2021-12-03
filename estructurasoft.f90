!factor de estructura de baxter
!
program baxter
implicit none
Integer, Parameter      :: nkmax = 2**12, nrmax=1000,nrmax2=500
Real(8), Parameter      ::  k0 = .01d0, R=1.d0
Real(8)			:: a, b, etha, rho, k, lamdak, ethavw, Ts, ri, sumablip, lm, x0, x1, xi ,integ, dr, fihs, dk, pi, dr1,ri1
Real(8), Dimension(nkmax) :: sdkArray
Real(8), Dimension(nrmax) :: us
Real(8)			:: sumagr
Real(8), Dimension(nrmax2) ::rd1,gdr
Complex(8), Parameter	:: I = Cmplx(0.d0,1.d0)
Integer		        :: nk, rdi, mu, ls,rdi1

pi=4.d0*atan(1.)

!!!!!!!!!!!11 integral de la blip

mu=6
Ts=0.5d0


! do ls=1,10

!!!!!!!!! limites de inttegracion :O

x0=0.d0		!	limites de la integral, integracion por trapecio
x1=1.d0

dr = (x1-x0)/nrmax*1d0
xi = x0
integ = ( f(x0)/2 ) *dr

do rdi=1, nrmax-1
 xi = xi + dr
 integ = integ + f(xi) * dr
 
 write(30,*) xi, f(xi)
end do

	integ = integ + (f(x1)/2 )*dr

lm=1.d0/((1.d0 - 3.d0*integ)**(1.d0/3.d0))				!lm=O_ss/O_hd

fihs=0.63d0 !+ (ls-1)*0.0185d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11


etha = fihs*(1.d0/lm)**3					!fraccion de volumen equivalente al sistema de esferas duras, esta es la que entra en la estructura

 print*,  lm, etha

write(100,*) 1.d0/lm,  Ts, etha

!stop

dk = .01d0

!!!!!!!!!!
ethavw=etha-(etha**2)/16.0d0	
				!!
lamdak = (ethavw/etha)**(1.0d0/3.0d0)		!!correccion verlet

rho = (6.0d0*ethavw)/(pi*R**3)


a=(1.d0 + 2.d0*ethavw)/((1.d0 - ethavw)**2)
b=-(3.d0/2.d0)*(R*ethavw)/((1.d0- ethavw)**2)

Open(Unit=1, File = "S(k)_fin.dat", Status = "Unknown")

dk=(dk/lm)


Do nk=1, nkmax
	
	k = (k0 + (nk-1)*dk)*lamdak

	sdkArray(nk)= 1.d0/(Q(k)*Q(-k))

	write(1,*) (k0 + (nk-1)*dk)*lm, sdkArray(nk)


enddo



contains

function f(x)
 real(8),intent(in) :: x
 real(8) :: f
!definimos la funcion a integrar


 f =(x**2)* exp(-(1.d0/Ts)*( (1.d0/x)**(2*mu) - 2.d0/(x)**mu + 1.d0))
 
 
end function f



complex(8) function Q(k)
real(8)	:: k

Q=1.d0 - ((pi*rho)/(k**3))*(2.d0*exp(I*k*R)*(b*k + a*(I + k*R)) - I*((2.d0*b*k*(-I + k*R)) + a*(2.d0 + (k**2)*(R**2))))

end function Q

end program baxter
