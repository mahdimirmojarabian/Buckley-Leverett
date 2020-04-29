from sympy import *
from matplotlib import pyplot as plt
from numpy import *

#parameters
muw=1.0e-3
muo=1.5e-3
kro0=1.
krw0=1.
no=2.
nw=2.
sor=0.2
swc=0.16
phi=0.25
k=1.0e-12

A=3.02700618 #area of the inlet face 
q=1e-4#rate of the inlet,m**3/day
t=0.35122*24.*3600.#inject time, second
u = q*t/A

L=100.#length of the model,m

pl=10.e6#pressure on the out boundary,Pa


def sws(sw):
    return (sw-swc)/(1-swc-sor)
def sos(sw):
    return (1-sw-sor)/(1-swc-sor)
def dsws(sw):
    return 1/(1-swc-sor)
def krw(sw):
    return krw0*sws(sw)**nw
def kro(sw):
    return kro0*(1-sws(sw))**no
    #return kro0*sos(sw)**2.*(1-(1-sos(sw))**no)
def dkrw(sw):
    return (dsws(sw)*nw*krw0*sws(sw)**(nw-1))
def dkro(sw):
    return (-dsws(sw)*no*kro0*(1-sws(sw))**(no-1))
def fw(sw):
    return ((krw(sw)/muw)/(krw(sw)/muw+kro(sw)/muo))
def dfw(sw):
    return ((dkrw(sw)/muw*kro(sw)/muo-krw(sw)/muw*dkro(sw)/muo)/(krw(sw)/muw+kro(sw)/muo)**2.)

# define initial conditions
sw0 = swc # the reservoir is filled with oil and connate water
# solve the nl equation to find the shock front saturation
def f_shock(sw):
    return (dfw(sw)-(fw(sw0)-fw(sw))/(sw0-sw))
swf=Symbol('swf')   
sw_shock = solve(f_shock(swf), swf)

#sw_shock has mutiple solutions and only one is physical
for i in range(len(sw_shock)):
    if sw_shock[i] > 0 and sw_shock[i] < 1:
        sw_shock_phys=sw_shock[i]

sw_inj = 1-sor  # not to make things too complicated
# water saturation between injection and shock
sw_rf = arange(sw_inj, sw_shock_phys, -0.001)
sw_all = arange(swc, 1-sor, 0.001) # whole range of possible sw
xt_all = u/phi*dfw(sw_all) # solution over a whole possible range of sw
xt_rf = u/phi*dfw(sw_rf)
xt_shock = xt_rf[-1]


def swx(x):#get sw of the location
    for sw in arange(0.,1.,0.0001):
        if x < xt_shock:
            if abs(x-u/phi*dfw(sw))<0.1:
                if sw <= sw_shock_phys:
                    return sw
        else:
            return swc

#pressure distribution
p=[]
x_l=[]
x_l.append(L)
p.append(pl)
n=100
for i in range(n-1,0,-1):
    x=L/n*i
    x_l.append(x)
    dx=L/n
    temp=dx*q/(k*kro(swx(x))/muo+k*krw(swx(x))/muw)/A
    p.append(temp+p[-1])

#plot
fig,ax=plt.subplots(2,2)

#relative permeability curve
ax[0][0].plot(sw_all, krw(sw_all),label='$k_{rw}$') 
ax[0][0].plot(sw_all, kro(sw_all),label='$k_{r0}$') 
ax[0][0].set_xlim(0,1)
ax[0][0].set_ylim(0,1)
ax[0][0].set_xlabel(r"$S_w$") 
ax[0][0].set_ylabel("Relative permeabilities") 
ax[0][0].legend()

#fw and dfw curve vs sw
line1=ax[1][0].plot(sw_all, fw(sw_all),label=r'$f_w$')
ax_t=ax[1][0].twinx()
line2=ax_t.plot(sw_all,dfw(sw_all),'g',label=r'$\frac{df_w}{dS_w}$')
ax[1][0].set_xlim(0,1)
ax[1][0].set_xlabel(r"$S_w$") 
ax[1][0].set_ylabel(r"$f_w$")
ax_t.set_ylabel(r'$\frac{df_w}{dS_w}$')
lns=line1+line2
labs=[l.get_label() for l in lns]
ax[1][0].legend(lns,labs,loc=0)

#sw distribution 
ax[0][1].plot(xt_all,sw_all,linestyle='--')
x1=[]
y1=[]
for i in range(50):
    x1.append(xt_shock+(max(xt_all)-xt_shock)/50.*i)
    y1.append(sw0)
x=hstack((xt_rf,array(x1)))
y=hstack((sw_rf,array(y1)))
ax[0][1].plot(x,y)
ax[0][1].set_ylim(0,1)
ax[0][1].set_xlim(0,100)
ax[0][1].set_xlabel("x/m") 
ax[0][1].set_ylabel(r"$S_w$") 
ax[0][1].minorticks_on()
ax[0][1].grid(which='minor',axis='both')

#pressure distribution
ax[1][1].plot(x_l,p)
ax[1][1].set_xlabel("x/m") 
ax[1][1].set_ylabel("Pressure/Pa") 
ax[1][1].minorticks_on()
ax[1][1].grid(which='minor',axis='both')

plt.show()

#save sw and p distribution to disk file
with open("{:.2f}-sw-day.txt".format(t/24./3600.),"w") as f:
    for i in range(len(x)):
        f.write("{}\t{}\n".format(x[i],y[i]))

with open("{:.2f}-p-day.txt".format(t/24./3600.),"w") as f:
    for i in range(len(x_l)):
        f.write("{}\t{}\n".format(x_l[i],p[i]))
