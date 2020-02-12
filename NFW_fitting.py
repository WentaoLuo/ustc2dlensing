#!/home/wtluo/anaconda/bin/python2.7

import emcee
import corner
import numpy as np 
import scipy.optimize as opt 
import matplotlib.pyplot as plt
from scipy import signal

pi      =np.pi
rho_crit=9.47e-27      # Critical density of the universe in unit of M_sun/(Mpc/h)^3
                      # at z=0, 3*H_0/(8*pi*G).
omega_m =0.23         # Density fraction of dark matter.
c_light =300000.0     # Speed of light in unit of km/s.
G_n     =6.754e-11    # Gravitational constant.
H0      =67.3         # Hubble constant at z=0.
ckm     =3.24078e-20  # 1 km equals ckm Mpc/h.
ckg     =5.0e-31      # 1 kg equals ckg M_sun.  rhoc    =rho_crit*ckg*10e+9/ckm/ckm/ckm
rhoc    =rho_crit*ckg*10e+9/ckm/ckm/ckm
rhom    =rhoc*omega_m

fac   =(c_light*c_light*ckg)/(4.*pi*G_n*ckm)

#---NFW ESD model----------------------------
def ESD(logM,con,f,rx,ry):
	nx      =int(np.sqrt(len(rx))) 
	Rr      =np.zeros((nx,nx))
	ESD     =np.zeros((nx,nx))
	xnew    =rx.reshape(nx,nx)
	ynew    =ry.reshape(nx,nx)

	rvir    =(10.0**logM*3.0/200./rhom/pi)**(1./3.)
        rs      =rvir/con
	ra      =rvir/con/np.sqrt(f)
	rb      =rvir*np.sqrt(f)/con
	theta   =np.arctan(ynew/xnew)
	sina    =np.sin(theta)
	cosa    =np.cos(theta)
	rs      =(ra*rb)/np.sqrt(ra*ra*sina*sina+rb*rb*cosa*cosa)
	delt    =(200.0/3.0)*(con**3)/(np.log(1.0+con)-con/(1.0+con))
	
	Rr      =np.sqrt(xnew*xnew/f+ynew*ynew*f)
	#print rs	
	for i in range(nx):
		for j in range(nx):
			Amp     =2.0*rs[i][j]*delt*rhoc
			#Amp     =2.0*rs*delt*rhoc
			if Rr[i][j]/rs[i][j]<1.0:
			#if Rr[i][j]/rs<1.0:
				xx       =Rr[i][j]/rs[i][j]
				#xx       =Rr[i][j]/rs
				f        =(1./(xx*xx-1.0))*\
					  (1.0-np.log((1.0+np.sqrt(1.0-xx*xx))/xx)/np.sqrt(1.0-xx*xx))
				g        =(2./(xx*xx))*\
					   (np.log(xx/2.)+np.log((1.0+np.sqrt(1.0-xx*xx))/xx)\
					   /np.sqrt(1.0-xx*xx))
				ESD[i][j]=10e-14*Amp*(g-f)
				#print xx
			elif Rr[i][j]/rs[i][j]==1.0:
			#elif Rr[i][j]/rs==1.0:
				xx       =Rr[i][j]/rs[i][j]
				#xx       =Rr[i][j]/rs
				f        =1./3.
				g        =2.+2.*np.log(0.5)
				ESD[i][j]=10e-14*Amp*(g-f)
			elif Rr[i][j]/rs[i][j]>1.0:
			#elif Rr[i][j]/rs>1.0:
				xx       =Rr[i][j]/rs[i][j]
				#xx       =Rr[i][j]/rs
				f        =(1./(xx*xx-1.0))*\
					  (1.0-np.arctan(np.sqrt(xx*xx-1.0))\
				   	  /np.sqrt(xx*xx-1.0))
				g        =(2./(xx*xx))*\
					  (np.log(xx/2.)+np.arctan(np.sqrt(xx*xx-1.0))\
					  /np.sqrt(xx*xx-1.0))
				ESD[i][j]=10e-14*Amp*(g-f)


	return ESD
#---END of NFW ESD model---------------------------------------------------

#---Likelyhood, Prior, Posterior-------------------------------------------
def lnprior(theta):
	logM,con,f=theta
	if 12.0<logM<15.5 and 1.<con<7.0 and 0.0<f<1.0:
		return 0.0
	return -np.inf
def lnlike(theta,rx,ry,g,er):
	logM,con,f=theta
	nx      =int(np.sqrt(len(rx)))
	model   =ESD(logM,con,f,rx,ry)
	nx      =int(np.sqrt(len(rx)))
	inv_s2  =(1.0/(er*er)).reshape(nx,nx)

	return -0.5*(np.sum((g-model)**2*inv_s2-np.log(inv_s2)))

def lnprob(theta,rx,ry,g,er):
	lp=lnprior(theta)
	if not np.isfinite(lp):
		return -np.inf
	return lp+lnlike(theta,rx,ry,g,er)
#---END of Likelyhood, Prior, Posterior-----------------------------------------

Sig_crit=fac*1215.2/882.9/332.3
zl,zs   =0.28,0.46
andis   =882.9

fname    ='ESDrot_400kpc_40_no200kpc'
#fname   ='ESDrot_theta_corrected'
#fname   ='test'
#fname   ='ESDnorot_400kpc_40'
#fname   ='testmcmc'
data     =np.loadtxt(fname,unpack=True)

rx       =data[0][:]
ry       =data[1][:]
g1       =data[2][:]
er       =data[4][:]/100.

nx       =int(np.sqrt(len(rx)))
g        =(g1).reshape(nx,nx)
inv_s2   =(1.0/(er*er)).reshape(nx,nx)

logM_ini=14.0
con_ini =5.6
f_ini   =0.6
sigma   =0.04
gsz_ini =1.0
model   =ESD(logM_ini,con_ini,f_ini,rx,ry)
#plt.imshow(g.T,interpolation='Nearest')
#plt.colorbar()
#plt.show()

chi2    =lambda *args:-2.0*lnlike(*args)
results =opt.minimize(chi2,[logM_ini,con_ini,f_ini],args=(rx,ry,g,er))

Mh,cc,ff=results["x"]
#print Mh,cc,ff
ndim,nwalkers=3,100
pos          =[results["x"]+1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler      =emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(rx,ry,g,er))
sampler.run_mcmc(pos,500)

burnin =100
samples=sampler.chain[:,burnin:,:].reshape((-1,ndim))

m_mc, c_mc, f_mc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],axis=0)))
print m_mc,c_mc,f_mc

from matplotlib.ticker import MaxNLocator

"""fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 9))
axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
axes[0].axhline(logM_ini, color="#888888", lw=2)
axes[0].set_ylabel("$logM$")

axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].axhline(con_ini, color="#888888", lw=2)
axes[1].set_ylabel("$c$")

axes[2].plot(sampler.chain[:, :, 2].T, color="k", alpha=0.4)
axes[2].yaxis.set_major_locator(MaxNLocator(5))
axes[2].axhline(f_ini, color="#888888", lw=2)
axes[2].set_ylabel("$f$")"""

fig1    =corner.corner(samples,lables=["logM","c","f","gsize"])
plt.show()

"""#figure()
plt.xlim([10,1500.0])
plt.xscale("log")
plt.ylim([1,1000.0])
plt.yscale("log")
plt.xlabel('r kpc/h')
plt.ylabel('ESD')
plt.errorbar(r0*1000.,g0,yerr=er0,marker='o',ms=8,color='black')
plt.plot(r0*1000,model0,'r-')
plt.show()"""

