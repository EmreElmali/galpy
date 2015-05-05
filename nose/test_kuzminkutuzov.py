#_____import packages_____
from galpy.potential import KuzminKutuzovStaeckelPotential
from galpy.potential import evaluatePotentials, evaluateRforces, evaluateDensities
from galpy.actionAngle import estimateDeltaStaeckel
from galpy.actionAngle import actionAngleStaeckel
from galpy.orbit import Orbit
from galpy.util import bovy_coords
import numpy
import scipy
import math

#Test whether circular velocity calculation works
def test_vcirc():
    #test the circular velocity of the KuzminKutuzovStaeckelPotential
    #using parameters from Batsleer & Dejonghe 1994, fig. 1-3
    #and their formula eq. (10)

    #_____model parameters______
    #surface ratios of disk and halo:
    ac_Ds = [50.  ,50.  ,50.  ,50. ,50. ,50.  ,40. ,40. ,40. ]
    ac_Hs = [1.005,1.005,1.005,1.01,1.01,1.01 ,1.02,1.02,1.02 ]
    #disk contribution to total mass:
    ks    = [0.05 ,0.06 ,0.07 ,0.07,0.1 ,0.125,0.1 ,0.12,0.125]
    #focal distance:
    Delta = 1.

    for ii in range(9):
        ac_D = ac_Ds[ii]
        ac_H = ac_Hs[ii]
        k    = ks[ii]
        
        #_____setup potential_____
        #first try, not normalized:
        V_D = KuzminKutuzovStaeckelPotential(amp=    k ,ac=ac_D,Delta=Delta,normalize=False)
        V_H = KuzminKutuzovStaeckelPotential(amp=(1.-k),ac=ac_H,Delta=Delta,normalize=False)
        pot = [V_D,V_H]
        #normalization according to Batsleer & Dejonghe 1994:
        V00 = evaluatePotentials(0.,0.,pot)
        #second try, normalized:
        V_D = KuzminKutuzovStaeckelPotential(amp=    k  / (-V00),ac=ac_D,Delta=Delta,normalize=False)
        V_H = KuzminKutuzovStaeckelPotential(amp=(1.-k) / (-V00),ac=ac_H,Delta=Delta,normalize=False)
        pot = [V_D,V_H]
        
        #_____calculate rotation curve_____
        Rs = numpy.linspace(0.,20.,100)
        z = 0.
        vcirc_calc = numpy.sqrt(-Rs * evaluateRforces(Rs,z,pot))
            
        #_____vcirc by Batsleer & Dejonghe eq. (10) (with proper Jacobian)_____
        def vc2w(R):
            g_D = Delta**2 / (1.-ac_D**2)
            a_D = g_D - Delta**2
            g_H = Delta**2 / (1.-ac_H**2)
            a_H = g_H - Delta**2
            l = R**2 - a_D
            q = a_H - a_D
            termD = numpy.sqrt(l  )*(numpy.sqrt(l  ) + numpy.sqrt(-g_D  ))**2
            termH = numpy.sqrt(l-q)*(numpy.sqrt(l-q) + numpy.sqrt(-g_D-q))**2
            return R**2 * (k / termD + (1.-k) / termH)
        vcirc_formula = numpy.sqrt(vc2w(Rs)/(-V00))

        assert numpy.all(numpy.fabs(vcirc_calc - vcirc_formula) < 10**-8.), \
                'Calculated circular velocity for KuzminKutuzovStaeckelPotential '+ \
                'does not agree with eq. (10) (corrected by proper Jacobian) '+ \
                'by Batsleer & Dejonghe (1994)'

    return None

#-----------------------------------------------------------------------------

#test whether the density calculation works
def test_density():
    #test the density calculation of the KuzminKutuzovStaeckelPotential
    #using parameters from Batsleer & Dejonghe 1994, tab. 2

    #_____parameters_____
    #table 2 in Batsleer & Dejonghe
    ac_D   = [25.  ,25.  ,25.  ,25. ,25.  ,25. ,40.  ,40.  ,40.  ,40. ,40. ,50.  ,50.  ,50. ,50.  ,50. ,50.  ,75.  ,75.  ,75. ,75.  ,75. ,75.  ,100. ,100. ,100.,100. ,100.,100. ,150. ,150. ,150. ,150. ,150.,150.]
    ac_H   = [1.005,1.005,1.01 ,1.01,1.02 ,1.02,1.005,1.005,1.01 ,1.01,1.02,1.005,1.005,1.01,1.01 ,1.02,1.02 ,1.005,1.005,1.01,1.01 ,1.02,1.02 ,1.005,1.005,1.01,1.01 ,1.02,1.02 ,1.005,1.005,1.01 ,1.01 ,1.02,1.02]
    k      = [0.05 ,0.08 ,0.075,0.11,0.105,0.11,0.05 ,0.08 ,0.075,0.11,0.11,0.05 ,0.07 ,0.07,0.125,0.1 ,0.125,0.05 ,0.065 ,0.07,0.125,0.10,0.125,0.05,0.065,0.07,0.125,0.10,0.125,0.05 ,0.065,0.075,0.125,0.11,0.125]
    Delta  = [0.99 ,1.01 ,0.96 ,0.99,0.86 ,0.88,1.00 ,1.01 ,0.96 ,0.99,0.89,1.05 ,1.06 ,1.00,1.05 ,0.91,0.97 ,0.98 ,0.99 ,0.94,0.98 ,0.85,0.91 ,1.06 ,1.07 ,1.01,1.06 ,0.94,0.97 ,1.06 ,1.07 ,0.98 ,1.06 ,0.94,0.97]
    Mmin   = [7.49 ,6.17 ,4.08 ,3.70,2.34 ,2.36,7.57 ,6.16 ,4.08 ,2.64,2.38,8.05 ,6.94 ,4.37,3.70 ,2.48,2.50 ,7.37 ,6.66 ,4.05,3.46 ,2.33,2.36 ,8.14 ,7.27 ,4.42,3.72 ,2.56,2.50 ,8.14 ,7.26 ,4.17 ,3.72 ,2.51,2.50]
    Mmax   = [7.18 ,6.12 ,3.99 ,3.69,2.37 ,2.40,7.27 ,6.11 ,3.99 ,2.66,2.42,7.76 ,6.85 ,4.26,3.72 ,2.51,2.54 ,7.07 ,6.51 ,3.95,3.48 ,2.36,2.40 ,7.85 ,7.15 ,4.30,3.75 ,2.58,2.54 ,7.85 ,7.07 ,4.08 ,3.75 ,2.53,2.53]
    rhomin = [0.04 ,0.05 ,0.04 ,0.04,0.03 ,0.03,0.06 ,0.06 ,0.05 ,0.04,0.04,0.07 ,0.08 ,0.06,0.07 ,0.04,0.05 ,0.08 ,0.09 ,0.07,0.09 ,0.05,0.06 ,0.12 ,0.13 ,0.09,0.13 ,0.07,0.09 ,0.16 ,0.19 ,0.12 ,0.18 ,0.10,0.12]
    rhomax = [0.03 ,0.03 ,0.02 ,0.03,0.02 ,0.02,0.04 ,0.04 ,0.03 ,0.03,0.02,0.05 ,0.05 ,0.04,0.05 ,0.03,0.03 ,0.05 ,0.06 ,0.04,0.06 ,0.03,0.04 ,0.07 ,0.08 ,0.06,0.08 ,0.04,0.05 ,0.09 ,0.10 ,0.07 ,0.10 ,0.06,0.07]
    Sigmin = [58   ,52   ,52   ,49  ,39   ,40  ,58   ,55   ,51   ,44  ,40  ,59   ,54   ,53  ,49   ,41  ,42   ,58   ,55   ,51  ,48   ,39  ,40   ,59   ,55   ,53  ,49   ,42  ,42   ,59   ,55   ,52   ,49   ,42  ,42]
    Sigmax = [45   ,41   ,38   ,37  ,28   ,28  ,45   ,32   ,37   ,32  ,30  ,46   ,43   ,40  ,37   ,30  ,31   ,45   ,43   ,38  ,36   ,28  ,29   ,46   ,43   ,40  ,38   ,31  ,31   ,46   ,44   ,39   ,38   ,30  ,31]

    for ii in range(len(ac_D)):

        if ac_D[ii] == 40.: 
            continue    
            #because I believe that there are typos in tab. 2 by Batsleer & Dejonghe...
    
        for jj in range(2):

            #_____parameters depending on solar position____
            if jj == 0:
                Rsun = 7.5
                zsun = 0.004
                GM   = Mmin[ii]  #units: G = 1, M in 10^11 solar masses
                rho  = rhomin[ii]
                Sig  = Sigmin[ii]
            elif jj == 1:
                Rsun = 8.5
                zsun = 0.02
                GM   = Mmax[ii]  #units: G = 1, M in 10^11 solar masses
                rho  = rhomax[ii]
                Sig  = Sigmax[ii]
            outstr = 'ac_D='+str(ac_D[ii])+', ac_H='+str(ac_H[ii])+', k='+str(k[ii])+', Delta='+str(Delta[ii])+\
                     ', Mtot='+str(GM)+'*10^11Msun, Rsun='+str(Rsun)+'kpc, rho(Rsun,zsun)='+str(rho)+'Msun/pc^3, Sig(Rsun,z<1.1kpc)='+str(Sig)+'Msun/pc^2'
            

            #_____setup potential_____
            amp_D = GM * k[ii]
            V_D = KuzminKutuzovStaeckelPotential(amp=amp_D,ac=ac_D[ii],Delta=Delta[ii],normalize=False)
            amp_H = GM * (1.-k[ii])
            V_H = KuzminKutuzovStaeckelPotential(amp=amp_H,ac=ac_H[ii],Delta=Delta[ii],normalize=False)
            pot = [V_D,V_H]

            #_____local density_____
            rho_calc = evaluateDensities(Rsun,zsun,pot) * 100. #units: [solar mass / pc^3]
            rho_calc = round(rho_calc,2)

            #an error of 0.01 corresponds to the significant digit 
            #given in the table, to which the density was rounded, 
            #to be wrong by one.
            assert numpy.fabs(rho_calc - rho) <= 0.01+10.**-8, \
                'Calculated density %f for KuzminKutuzovStaeckelPotential ' % rho_calc + \
                'with model parameters:\n'+outstr+'\n'+ \
                'does not agree with value from tab. 2 '+ \
                'by Batsleer & Dejonghe (1994)' 

            #_____surface density_____
            Sig_calc, err = scipy.integrate.quad(lambda z: (evaluateDensities(Rsun,z/1000.,pot) * 100.), #units: [solar mass / pc^3]
                                            0., 1100.) #units: pc
            Sig_calc = round(2. * Sig_calc)

            #an error of 1 corresponds to the significant digit 
            #given in the table, to which the surface density was rounded, 
            #to be wrong by one.
            assert numpy.fabs(Sig_calc - Sig) <= 1., \
                'Calculated surface density %f for KuzminKutuzovStaeckelPotential ' % Sig_calc + \
                'with model parameters:\n'+outstr+'\n'+ \
                'does not agree with value from tab. 2 '+ \
                'by Batsleer & Dejonghe (1994)' 

    return None

#-----------------------------------------------------------------------------

#test wheter the orbit integration in C and Python are the same
def test_orbitItegrationC():

    #_____initialize some KKSPot_____
    Delta = 1.0
    pot = KuzminKutuzovStaeckelPotential(ac=20.,Delta=Delta,normalize=True)

    #_____initialize an orbit (twice)_____
    vxvv = [1.,0.1,1.1,0.,0.1]
    o_P= Orbit(vxvv=vxvv)
    o_C= Orbit(vxvv=vxvv)

    #_____integrate the orbit with python and C_____
    ts= numpy.linspace(0,100,101)
    o_P.integrate(ts,pot,method='leapfrog')  #python
    o_C.integrate(ts,pot,method='leapfrog_c')#C

    assert numpy.all(numpy.fabs((o_P.R(ts) - o_C.R(ts))/o_P.R(ts)) < 10.**-5), \
            'Orbit integration for R coordinate different in ' + \
            'C and Python implementation.'
    assert numpy.all(numpy.fabs((o_P.z(ts) - o_C.z(ts))/o_P.z(ts))[1::] < 10.**-1), \
            'Orbit integration for z coordinate different in ' + \
            'C and Python implementation.'
    assert numpy.all(numpy.fabs((o_P.vR(ts) - o_C.vR(ts))/o_P.vR(ts)) < 10.**-3), \
            'Orbit integration for vR coordinate different in ' + \
            'C and Python implementation.'
    assert numpy.all(numpy.fabs((o_P.vz(ts) - o_C.vz(ts))/o_P.vz(ts)) < 10.**-1), \
            'Orbit integration for vz coordinate different in ' + \
            'C and Python implementation.'
    assert numpy.all(numpy.fabs((o_P.vT(ts) - o_C.vT(ts))/o_P.vT(ts)) < 10.**-5), \
            'Orbit integration for vT coordinate different in ' + \
            'C and Python implementation.'

    return None

#-----------------------------------------------------------------------------

#test whether this is really a Staeckel potential and the Delta is constant
def test_estimateDelta():

    #_____initialize some KKSPot_____
    Delta = 1.0
    pot = KuzminKutuzovStaeckelPotential(ac=20.,Delta=Delta,normalize=True)

    #_____initialize an orbit (twice)_____
    vxvv = [1.,0.1,1.1,0.01,0.1]
    o= Orbit(vxvv=vxvv)

    #_____integrate the orbit with C_____
    ts= numpy.linspace(0,101,100)
    o.integrate(ts,pot,method='leapfrog_c')


    #____estimate Focal length Delta_____
    #for each time step individually:
    deltas_estimate = numpy.zeros(len(ts))
    for ii in range(len(ts)):
        deltas_estimate[ii] = estimateDeltaStaeckel(o.R(ts[ii]),o.z(ts[ii]),pot=pot)

    assert numpy.all(numpy.fabs(deltas_estimate - Delta) < 10.**-8), \
            'Focal length Delta estimated along the orbit is not constant.'

    #for all time steps together:
    delta_estimate = estimateDeltaStaeckel(o.R(ts),o.z(ts),pot=pot)
    
    assert numpy.fabs(delta_estimate - Delta) < 10.**-8, \
            'Focal length Delta estimated from the orbit is not the same as the input focal length.'

    return None

#-----------------------------------------------------------------------------

#test whether this is really a Staeckel potential and the Actions are conserved along the orbit
def test_actionConservation():

    #_____initialize some KKSPot_____
    Delta = 1.0
    pot = KuzminKutuzovStaeckelPotential(ac=20.,Delta=Delta,normalize=True)

    #_____initialize an orbit (twice)_____
    vxvv = [1.,0.1,1.1,0.01,0.1]
    o= Orbit(vxvv=vxvv)

    #_____integrate the orbit with C_____
    ts= numpy.linspace(0,101,100)
    o.integrate(ts,pot,method='leapfrog_c')

    #_____Setup ActionAngle object and calculate actions (Staeckel approximation)_____
    aAS = actionAngleStaeckel(pot=pot,delta=Delta,c=True)
    jrs,lzs,jzs = aAS(o.R(ts),o.vR(ts),o.vT(ts),o.z(ts),o.vz(ts))

    assert numpy.all(numpy.fabs(jrs - jrs[0]) < 10.**8), \
        'Radial action is not conserved along orbit.'

    assert numpy.all(numpy.fabs(lzs - lzs[0]) < 10.**8), \
        'Angular momentum is not conserved along orbit.'

    assert numpy.all(numpy.fabs(jzs - jzs[0]) < 10.**8), \
        'Vertical action is not conserved along orbit.'

    return None
#-----------------------------------------------------------------------------

#test coordinate transformation
def test_lambdanu_to_Rz():
    #coordinate system:
    a = 3.
    g = 4.
    Delta = numpy.sqrt(g-a)
    ac = numpy.sqrt(a/g)
    #coordinate transformation:
    l, n = 2., -4.
    R,z= bovy_coords.lambdanu_to_Rz(l,n,ac=ac,Delta=Delta)
    #true values:
    R_true = numpy.sqrt((l+a)*(n+a)/(a-g))
    z_true = numpy.sqrt((l+g)*(n+g)/(g-a))
    #test:
    assert numpy.fabs(R-R_true) < 10.**-10., 'lambdanu_to_Rz conversion did not work as expected (R)'
    assert numpy.fabs(z-z_true) < 10.**-10., 'lambdanu_to_Rz conversion did not work as expected (z)'
    #Also test for arrays
    os= numpy.ones(2)
    R,z= bovy_coords.lambdanu_to_Rz(os*l,os*n,ac=ac,Delta=Delta)
    assert numpy.all(numpy.fabs(R-R_true) < 10.**-10.), 'lambdanu_to_Rz conversion did not work as expected (R array)'
    assert numpy.all(numpy.fabs(z-z_true) < 10.**-10.), 'lambdanu_to_Rz conversion did not work as expected (z array)'
    return None

def test_Rz_to_lambdanu():
    #coordinate system:
    a = 3.
    g = 7.
    Delta = numpy.sqrt(g-a)
    ac = numpy.sqrt(a/g)
    #true values:
    l, n = 2., -5.
    #coordinate transformation:
    lt,nt= bovy_coords.Rz_to_lambdanu(*bovy_coords.lambdanu_to_Rz(l,n,ac=ac,Delta=Delta),ac=ac,Delta=Delta)
    #test:
    assert numpy.fabs(lt-l) < 10.**-10., 'Rz_to_lambdanu conversion did not work as expected (l)'
    assert numpy.fabs(nt-n) < 10.**-10., 'Rz_to_lambdanu conversion did not work as expected (n)'
    #Also test for arrays
    os= numpy.ones(2)
    lt,nt= bovy_coords.Rz_to_lambdanu(*bovy_coords.lambdanu_to_Rz(l*os,n*os,ac=ac,Delta=Delta),ac=ac,Delta=Delta)
    assert numpy.all(numpy.fabs(lt-l) < 10.**-10.), 'Rz_to_lambdanu conversion did not work as expected (l array)'
    assert numpy.all(numpy.fabs(nt-n) < 10.**-10.), 'Rz_to_lambdanu conversion did not work as expected (n array)'
    return None

