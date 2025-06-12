import numpy as np
from astropy.io import ascii
import scipy.io
import scipy.integrate
import scipy.special as funcs
import numpy as np
import scipy.integrate as integrate
import radxglobals as cnst



def prep_pmesh(z):
# z is an irregular grid and must be at cell boundaries for pcolormesh (therefore make an array that is ndep + 1 dimensions.)
    ndep = len(z)
    midz = (z[1:len(z)] + z[0:len(z)-1])/2.
    newz = np.insert(midz, 0, z[0] + (z[0]-midz[0]))
    ndep2=len(newz)
    z_bdry = np.append(newz, z[ndep-1] + (z[ndep-1]-midz[ndep-2]))
    return z_bdry



def findind(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def monoEsy(E_MeV = 85.0, B_gauss = 15000., pitch_angle = 0.0):
    ''' ret_wl, ret_spec, integral1, integral2, rL_cm, lam_c = ntspec.monoEsy(E_MeV = 85.0, B_gauss = 15000., pitch_angle = 0.deg)   
    E_MeV is Kinteic Energy only! 15000 T agrees with NIST radius of curvature of 0.83 m for E=380 MeV
    Calculates synchrotron spectrum in erg/s/cm2/Ang for mono-energetic electron spiraling in magnetic field.  
    Returns, wl(angstrom), spectral power irradiance(erg/s/Ang) per electron, integral of spectral power irradiance, standard rybicki and lightman formula for power (if they don't agree, need to modify wl_array below), radius of orbit in cm, and the characteristic wavelength (angstroms), sometimes also called the critical photon energy 
 '''
    Energy = E_MeV # MeV
    rest_mass = 0.511
    theta_array1 = np.linspace(0,0.511 / Energy, 20) * 180. / cnst.PI
    theta_array2 = np.linspace(0.511 / Energy *1.2, 0.511 / Energy * 5, 5) * 180. / cnst.PI
    theta_array = np.concatenate([theta_array1,theta_array2])
    wl_array=10**np.arange(0, 4.5, 0.01)
    Plam = calc_Plam_sync(wl_array,theta_array, Energy_MeV=Energy, Bfield_G = B_gauss)
    integ = np.trapz(Plam, x=wl_array)

   # mass = 9.109e-28 * rest_mass / 0.511 / 1000.0
    
    relectron = 2.8179409e-13 # cm
    gamma = E_MeV / rest_mass + 1.0
    beta = np.cos(pitch_angle/180. * np.pi) * np.sqrt(1. - 1./gamma**2)
    rho = gamma * beta * cnst.MEL * cnst.CC**2 / ( B_gauss * cnst.EL)  # need an extra factor of c for cgs.
    integ_Wied = (relectron/100.0) * (cnst.MEL /1000.0) * (cnst.CC/100.0)**3 * gamma**4 / (rho/100.0)**2 * 2./3. * 1e7 # convert to erg/s.  The two integrals match!
    phot_crit = 1.5 * 1.054e-34 * gamma**3 * cnst.CC / 100.0 / (rho / 100.0) * 1e7  # Wiedemann 2003 eq 9.78.
 #   print(phot_crit/1e7*6.24e15 / (0.66503 * (E_MeV/1000.)**2 * (B_gauss/1e4)))
    lam_max_ang = cnst.HPL * cnst.CC / phot_crit * 1e8  # this is lam_c in NIST SURF calculator!

    indpeak = np.argmax(Plam)
    lam_peak_ang = wl_array[indpeak]
   # in per angstrom, plot,x,exp(-1./x) * x^(-2.5) where x = lam/lamc and x goes from 0 to 1.
    return wl_array, Plam, integ, integ_Wied, rho, lam_max_ang, lam_peak_ang

def calc_Plam_sync(wl_array, theta_array, Energy_MeV=180.0, Bfield_G=15000.0):
    Plam_rad = np.zeros((len(wl_array),len(theta_array)))
    Plam = np.zeros(len(wl_array))
    ww=0
    while ww < len(wl_array):
        rr = 0
        while rr < len(theta_array):
         #   if Bfield_G < 0:
         #       Plam_rad[ww, rr] = sync_spec_rad(wl_array[ww], theta_array[rr], E_user = Energy_MeV)
            if Bfield_G > 0:
                Plam_rad[ww, rr] = sync_spec_rad(wl_array[ww], theta_array[rr], E_user = Energy_MeV, B_user=Bfield_G)
            rr=rr+1
        Plam[ww] = 2.0 * scipy.integrate.trapz(Plam_rad[ww,:],x=theta_array /180. * cnst.PI)
         
        ww=ww+1
    return Plam
    
def sync_spec_rad(user_lam, phi_degree, E_user = 180.0, B_user = -100, return_option = 'P'):
    # user_lam in Angstroms, E_user is kinetic energy in MeV, phi_degree in degrees, B_user in Gauss
    
    Ein = E_user
    rest_mass = 0.511
    mass = 9.109e-28 * rest_mass / 0.511 / 1000.0
    q = 1.602e-19  # 4.803e-10 esu (1 C = 3e9 esu)
    phi = phi_degree / 180.0 * cnst.PI
    gamma = E_user / rest_mass + 1.0
    ### Need to add pitch angle to vperp!
    if B_user > 0:
        B_Tesla = B_user / 1e4
        vperp = ((1.0 - 1./ (Ein/(rest_mass) + 1.0)))**0.5 * 2.998e8
        Rorbit = vperp / (q * B_Tesla) * mass * 100.0 * gamma  # cm
        
    if B_user < 0:
        Rorbit = 83.4 * (Ein/180.)**0.5 # cm  # scaled to NIST orbit.
        vperp = ((1.0 - 1./ (Ein/(rest_mass) + 1.0)))**0.5 * 2.998e8
        B_Gauss = vperp / (q * Rorbit/100) * mass * 1e4 # gyroradius in SI units (meters)
        
  #  print('The rL in cm is', Rorbit)

    #phi = rest_mass / Ein
   # print(phi * 180./cnst.PI, 'degrees')
    lamin = user_lam / 1e8  # Ang to cm.
    
    x = Ein / rest_mass * phi
    xi = 2.0 * cnst.PI * Rorbit / (3.0 * lamin) * (Ein/rest_mass)**(-3.0) * (1+x**2)**(3./2.)
    f_o = cnst.CC / Rorbit
    order = 2./3.
    K_23 = funcs.kv(order, xi) # modified bessel functions of the second kind.
    order = 1./3.
    K_13 = funcs.kv(order, xi)
    P_lam = 8.0/3.0 * cnst.PI * cnst.EL**2 * cnst.CC**2 / ( f_o * lamin**4) * (rest_mass / Ein)**4 * \
        (1.0 + x**2)**2 * (K_23**2 + K_13**2 * (x**2 / (1.0 + x**2))) * 1e-8  # This last factor converts 
    # from /cm to /Ang, so I_lam is in units of erg/s/radian/Ang.

    # return values all in cgs!
    if return_option == 'P':
        ret_val = P_lam
    if return_option == 'B':
        ret_val = B_Gauss
    if return_option == 'L':
        ret_val = Rorbit
        
    return ret_val


def calc_snu_alpha(frq1, frq2, iHz1, iHz2):
    alpha = (np.log10(iHz2) - np.log10(iHz1))/(np.log10(frq2) - np.log10(frq1)) #  positive alpha increasing with frequency (lower frequency is frq1)
    delta_thin  = (1.22 - alpha)/0.9  # power-law index assuming optically thin emission

    return alpha, delta_thin


def gyrosynchPL(Hz_grid, Bgauss=50.0, in_theta=45.,delta=3.0, N10keV=1e10, L = 3.4e8, N20keV = -9):
    # optically thin case.
    # N10keV is calcualted in NT bremss routine. See dulk 1985.

    if N20keV > 0:
        N10keV = N20keV * (20.0 / 10.0)**(delta)
    
    # nu_B, wave_grid, jlam, Nphot, nu_peak, jHz, tau, iHz = nt.gyrosynchPL(np.array([34e9, 17e9]), delta = 4.7, Bgauss=800.0)
    theta = in_theta / 180.0 * np.pi
    nu_B = 2.8e6 * Bgauss
    lam_B_Ang = cnst.CCANG / nu_B
    wave_grid = cnst.CC / Hz_grid * 1e8
    # dulk's emissivity eta_nu is erg/s/cm3/Hz/sr.
    jlam = Bgauss * 3.3e-24 * 10**(-0.52 * delta) * (np.sin(theta))**(-0.43 + 0.65*delta) * ( Hz_grid / nu_B )**(1.22-0.90*delta) * cnst.CCANG / (wave_grid)**2 * N10keV
    jHz = Bgauss * 3.3e-24 * 10**(-0.52 * delta) * (np.sin(theta))**(-0.43 + 0.65*delta) * ( Hz_grid / nu_B )**(1.22-0.90*delta) * N10keV
    # jlam in erg/s/cm3/sr/Ang, jHz is  eta_nu is erg/s/cm3/Hz/sr
    Nphot = np.trapz(jlam / (cnst.CC * cnst.HPL/(wave_grid/1e8)), x = wave_grid)

   # L = (10. + 3. + 1.)/3. * 735e5  # average of dims in White 2011
    nu_peak = 2.72 * 10**(3.0 + 0.27*delta) * (np.sin(theta))**(0.41 + 0.03*delta) * (N10keV*L)**(0.32 - 0.03*delta) * Bgauss**(0.68+0.03*delta)

    kappa = N10keV/Bgauss * 1.4e-9*10**(-0.22*delta)*(np.sin(theta))**(-0.09 + 0.72*delta) * (Hz_grid / nu_B)**(-1.30-0.98*delta)
    Teff = 2.2e9 * 10**(-0.31 * delta) * (np.sin(theta))**(-0.36 - 0.06*delta) * (Hz_grid/nu_B)**(0.50 + 0.085 * delta)
    
    tau = kappa * L
    Src = jHz/kappa
    iHz =  Src * (1.0-np.exp(-tau))
   # ilam = Io * np.exp(-tau_line) + Src * (1.0-np.exp(-tau_line)) # erg/s/cm2/Ang/sr emergent from a simple slab.
    return nu_B, wave_grid, jlam, Nphot, nu_peak, jHz, tau, iHz, Teff

def gyrorad(Bgauss, KE=37.0, mass=9.1093897E-28,utheta=45.0,Z=1.0):
    ''' your_dict = gyrorad(Bgauss, KE=37, mass=9.109e-28, utheta=45.0,Z=1.0) 
    KE  is in keV.
'''
    #mp =1.6726231E-24
    rest_nrg = mass * cnst.CC**2 * 1./cnst.KEV2ERG # keV
    gamma = KE / rest_nrg + 1.0 # KE = (gamma-1)*mc**2
    vel = cnst.CC * np.sqrt( (gamma**2 - 1.0) /gamma**2 ) # cm/s
    vperp  = vel * np.cos(utheta/180.0 * np.pi)
    rL_cm = mass * vperp *  cnst.CC / (Z*cnst.EL * Bgauss) * gamma
    freqB = Z * cnst.EL * Bgauss / (mass * cnst.CC) / 2.0 / np.pi   #= 2.8e6  * Bgauss for electrons...  # revolutions / sec.
    freqB_rel =  vperp / rL_cm / 2.0 / np.pi
    freq = freqB_rel * gamma**3
    vz = vel * np.sin(utheta/180. * np.pi)
    tflight = 1e9 / vz
    revs10Mm = tflight * freqB_rel
    
    udict = {'beta':vel/cnst.CC, 'rL_cm':rL_cm, 'gamma':gamma, 'revGHz':freqB/1e9, 'tflight_s_10Mm':tflight,  'revs10Mm':revs10Mm, 'harmonic_GHz':freq/1e9, 'freqB_rel':freqB_rel}
    return udict


def rftab(ftabname):
    f = open(ftabname, 'r')

    E = []
    F = []
    for line in f:
        linef = line
        linef = linef.strip()
        linef = linef.split()
        if linef[0] == '*Injected':
            if linef[4] == 'boundaries':
                linee = f.readline().strip().split()
                while linee[0] != '*Time':
                   # print('found it ******************')
                    E = np.append(E, np.array(linee))
                    linee = f.readline().strip().split()
                t = f.readline()
                t = f.readline()
                lineflx = f.readline().strip().split()
                while lineflx[0] != '*Time':
                    F = np.append(F, np.array(lineflx))
                    lineflx = f.readline().strip().split()
    EE = np.array(E,dtype='float64')
    FF = np.array(F,dtype='float64')

    return EE, FF

def synchPL(Hz_grid, Bgauss=50.0, in_theta=45.,delta=3.0, N10MeV=1e10, Eo=10000.):
    # optically thin case.
    # N10keV is calcualted in NT bremss routine. See dulk 1985.
    # Need to check how high energies of electrons go in MeV.
    Ecut = Eo
    theta = in_theta / 180.0 * np.pi
    nu_B = 2.8e6 * Bgauss
    lam_B_Ang = cnst.CCANG / nu_B
    wave_grid = cnst.CC / Hz_grid * 1e8
    # dulk's emissivity eta_nu is erg/s/cm3/Hz/sr.
    jlam = Bgauss * 8.6e-24*(delta-1.0) * np.sin(theta) * (0.175 / np.sin(theta) * (Ecut/ 1000.)**(-2.) * Hz_grid/nu_B)**(-1.0*(delta-1)/2.)  * cnst.CCANG / (wave_grid)**2 * N10MeV
    # jlam in erg/s/cm3/sr/Ang
    Nphot = np.trapz(jlam / (cnst.CC * cnst.HPL/(wave_grid/1e8)), x = wave_grid)

    jHz =  Bgauss * 8.6e-24*(delta-1.0) * np.sin(theta) * (0.175 / np.sin(theta) * (Ecut/ 1000.)**(-2.) * Hz_grid/nu_B)**(-1.0*(delta-1)/2.) * N10MeV
    L = (10. + 3. + 1.)/3. * 735e5  # average of dims in White 2011
    nu_peak = 3.2e7 * np.sin(theta) * (Ecut/1000.)**( (2.0 * delta -2.)/(delta+4.0)) * (8.7e-12 * (delta-1.0)/np.sin(theta) * N10MeV * L)**(2./(delta+4.0)) * Bgauss**((delta+2.0)/(delta+4.0))

    
    kappa = N10MeV/Bgauss * 8.7e-12 * (delta-1.)/np.sin(theta) * (Ecut/1000.0)**(delta-1.0)*(8.7e-2/np.sin(theta) * Hz_grid/nu_B)**((-1.0)*(delta+4.0)/2.)
    tau = kappa * L
    Src = jHz/kappa
    iHz =  Src * (1.0-np.exp(-tau))
    return nu_B, wave_grid, jlam, Nphot, nu_peak, jHz,tau, iHz


def WWBW(scale=1.0, column = 2):
    dat = ascii.read('/home/adamkowalski/Desktop/Projects/langmuirbeam/Kontar2012_langmuir_HighRes')

    e0 = dat['x']

    if column == 2:
        blue = dat['Curve2(blue)']
    else:
        blue = dat['Curve1(solidblack)']
        
    print('Integrated flux = {0:.3e}'.format(np.trapz(10**blue, x=10**e0)))

    emesh = prep_pmesh(e0)

    print('Integrated flux after scalling = {0:.3e}'.format(np.trapz(scale*10**blue, x=10**e0)))
    print('Integrated energy flux after scalling = {0:.3e}'.format(np.trapz(scale*10**blue * 10**e0 * 1.602e-9, x=10**e0)))

    scaled_flux = scale * np.array(10**blue,dtype=float)

    integrated_flux = np.array(np.trapz(scale*10**blue, x=10**e0),dtype=float)

    return 10**e0, scaled_flux
  

def ntbremss(Ecut=37.0,Emax=50000.,delta=3.0,nEE=500.,Eflux=1e13, n_amb = 1e15, pasig=-99., Bg = 600.0, calc_cttm = 0, radio_L = 3.4e8, m_particle= 511.0, pasig_method='theta', read_ftab = '?', wwbw=0.0, delta_2 = 5.0):
    #c: just like in RADYN
    #c: delta is delta input in radyn so delta in number electrons / s /cm2/ keV
    #c: isotropic in forward hemisphere.
    ''' your_dict = ntspec.ntbremss(Ecut=37.0,delta=3.0,Eflux=1e13, n_amb = 1e15, pasig=-99.,Emax=5e4, nEE=500) '''
    # to calculate CTTM spectrum, set calc_cttm=1;  the dictionary variables you'll want are eps_arr (keV) and I_keV_thick (photons/s/cm2/keV)
    
    if pasig < 0:
        ibeam = 'iso'

    Emin = Ecut
    keV2erg = 1.602E-9
    erg2keV = 1./keV2erg
    gs=0.5
    #Emax=300000. # 300 MeV; in radyn, max is 50 MeV.
    kboltz=1.3806e-16
    ergperev = 1.60217657e-12
    nEbc = round(nEE/3)
    Emn = (5.0 * 3. * kboltz * 1e4 /ergperev *1e-3)**gs #c:  3 k T where T =1d4
    Emx = (Emin*0.998)**gs
    E = np.zeros(int(nEE))
    
   
    for j in range(int(nEbc)):
        E[j] = ((float(j))/(nEbc-1.)*(Emx - Emn) + Emn)**(1./gs)

   
    for j in range(int(nEE-nEbc)):
        E[j+nEbc]=np.exp((float(j))/(nEE-nEbc-1.)*(np.log(Emax)-np.log(Emin))+ np.log(Emin))

    #c: E is electron dist energy grid in [keV]   
    good = np.all([E >= Emin],axis=0)  #c: forcing a low-energy cutoff becuase I don't have thin target
    E_keV = E[good]
  #c:  test[1:-1]-test[0:-2]
    
    #c: delta is a power law in F_N = number flux [el / cm2 / s / keV], capital non-curly F in Holmanetal 2011 SSSRev.
    F_n = E_keV**(-1.0 * delta)
    A = np.trapz(F_n * E_keV * keV2erg, x=E_keV) / Eflux  #c: A * F_N * E * 1.602e-9 * dE is ergs/cm2/s/keV * dkeV
    #c: np.tsum(E_keV, F_N * E_keV) * A should give number / s /cm2
    #c:  I = A * F_erg / pi is the isotropic pitch angle distribution function.
    Eflux10 = Eflux * (10./Ecut)**(2.0-delta)
    g10 = np.all([E >=10.],axis=0)
    E10 = E[g10]
    F_n10 = E10**(-1.0 * delta)
    A10 = np.trapz(F_n10 * E10 * keV2erg, x=E10) / Eflux10  #c: A * F_N * E * 1.602e-9 * dE is ergs/cm2/s/keV * dkeV
    F_N10 = F_n10 / A10
    F10 = np.trapz(F_N10, x=E10)
    gamma10 = 1.0 + E10 / 511.0
    vel_E10 = cnst.CC * np.sqrt( (gamma10**2 - 1.0) /gamma10**2 )
    NN_E_iso10 = (F_N10 / np.pi) / vel_E10 * 2.0 * np.pi #c: this is number / cm3 / keV, similar to Rutten eq. 2.9 except integrating
    NN_iso10 = np.trapz(NN_E_iso10, x=E10)  # number el/cm3 for radio gyrosynch.  Need to calculate using 10keV cutoff.

    F_N = F_n * 1./A
    F = np.trapz(F_N, x=E_keV) # = NTel/cm2/s flux.
   # print(F, F_N)
    
    gamma = 1.0 + E_keV / m_particle
    vel_E = cnst.CC * np.sqrt( (gamma**2 - 1.0) /gamma**2 )
    NN_E_iso = (F_N / np.pi) / vel_E * 2.0 * np.pi #c: this is number / cm3 / keV, similar to Rutten eq. 2.9 except integrating
    # over only one hemisphere.

    if pasig >= 0.0:
        if read_ftab == 'other':
           # ftab_dat = ascii.read(read_ftab)
            E_keV_coarseedge, F_N_coarse = rftab(read_ftab)
            E_keV_coarse = E_keV_coarseedge[1:,]*0.5 + E_keV_coarseedge[0:-1]*0.5
            #E_keV_coarse = 10**np.array(ftab_dat['log10EkeV'])
            good_int = np.all([E_keV <= np.max(E_keV_coarse), E_keV >= np.min(E_keV_coarse)],axis=0)
            #F_N_coarse = 10**np.array(ftab_dat['log10El_scm2keV'])
            F_N = 10**np.interp(np.log10(E_keV[good_int]), np.log10(E_keV_coarse), np.log10(F_N_coarse))
            vel_E = vel_E[good_int]
            E_keV = E_keV[good_int]
            NN_E_iso = (F_N / np.pi) / vel_E * 2.0 * np.pi
        
        NRG_ARR, DNRG_ARR, FLUX_DIST, MU_ARR, DMU_ARR, THETA_ARR, NUM_DIST, INTENSITY_DIST, TOTNUMFLUX, FNUMBER, BEAMDENSITY = power_law(E_keV, vel_E, F_N * E_keV * keV2erg, NN_E_iso, PASIG=pasig, PASIG_METHOD=pasig_method) # energy(kev), relativistic velocity (cm/s), spectral energy flux (erg/s/cm2/keV), number density (number / cm3 / keV)
        return NRG_ARR, DNRG_ARR, FLUX_DIST, MU_ARR, DMU_ARR, THETA_ARR, NUM_DIST, INTENSITY_DIST, TOTNUMFLUX, FNUMBER, BEAMDENSITY


    # subscript E is for electron keV, epsilon is for photon keV
    NN_iso = np.trapz(NN_E_iso, x=E_keV)  # number el/cm3 for radio gyrosynch.  Need to calculate using 10keV cutoff:  NN_iso10
    NN_iso_E = np.zeros_like(E_keV)
    for ele in range(len(E_keV)-1):
        NN_iso_E[ele] = np.trapz(NN_E_iso[ele:,] ,x=E_keV[ele:,])
        
    
   # j_eps = np.tsum(E_keV, n_amb * Q_E * vel_E * NN_E, eps, Emax) / 4.0 / np.pi #c:  cm-3 * cm2/photkeV * cm/s (vel) * el/cm3/keV (NN) *  dkeV (dE_keV) / 4pi, emitted equally in all directions though beam was in one direction only.
    #c: this equation is in Holman et al. 2011 first paragerph of section 2.
    #  nu_eps is the photon_yield in Holman et al 2011 which is n_amb * Q_E * vel = photons(eps) / s / photkeV / electron(E) integrated over time


    j_eps = np.zeros(len(E))
    I_thick = np.zeros(len(E))
    F_NN = np.zeros(len(E))
    F_NN[good] = F_N

    if read_ftab == 'ftab.dat':
        EK12edge, F_NK12 = rftab(read_ftab)
      #  EK12, F_NK12 = WWBW(scale=wwbw)
        EK12 = 10**( np.log10(EK12edge[1:,]) * 0.5 + np.log10(EK12edge[0:-1]) * 0.5)
        F_N = 10**np.interp( np.log10(E), np.log10(EK12), np.log10(F_NK12))
        bad=np.any([E < 10.0, E > Emax],axis=0)
        F_N[bad] = 1e-30
        print('Reading in WWBW from FTAB.DAT')
       # if Emax > 430:
        #    EmaxK12 = EK12[-1]
       #     FmaxK12 = F_NK12[-1]
            #Flim = -1.0 * 10**(delta_2 * (np.log10(Emax) - np.log10(EmaxK12)) + np.log10(FmaxK12))
        #    good2 = (E > EmaxK12)
       #     F_N[good2] = FmaxK12 * (E[good2]/EmaxK12)**(delta_2 * (-1.0))
        F_NN = F_N
    

    
    if wwbw > 0.0:
        EK12, F_NK12 = WWBW(scale=wwbw)
        F_N = 10**np.interp( np.log10(E), np.log10(EK12), np.log10(F_NK12))
        bad=np.any([E < 10.0, E > Emax],axis=0)
        F_N[bad] = 1e-30
        print('Reading in WWBW from Kontar et al. 2012')
        if Emax > 430:
            EmaxK12 = EK12[-1]
            FmaxK12 = F_NK12[-1]
            #Flim = -1.0 * 10**(delta_2 * (np.log10(Emax) - np.log10(EmaxK12)) + np.log10(FmaxK12))
            good2 = (E > EmaxK12)
            F_N[good2] = FmaxK12 * (E[good2]/EmaxK12)**(delta_2 * (-1.0))
        F_NN = F_N
    
    Q_E = np.zeros((len(E),len(E)))
    Integ1 = np.zeros(len(E))
    eps_arr = E  # photon energy grid same as electron energy grid to << Ecut  [keV]
    for ee in range(len(eps_arr)):
        # This is non-relativistic formula cross section in cm2/photkev (see Holman et al. 2011):
        realE = np.all([E > eps_arr[ee]],axis=0)
        Q_E[ee,realE] = 1.4 * 7.90e-25 / eps_arr[ee] / E[realE] * np.log( (1.0 + np.sqrt(1.0 - eps_arr[ee]/E[realE])) / (1.0 - np.sqrt(1.0 - eps_arr[ee]/E[realE])) )  # cross section of production of photon energy epsilon as a function of electron energy greater than it.
        j_eps[ee] = np.trapz( n_amb * Q_E[ee, ee:,] * F_NN[ee:,], x = E[ee:,]) / 4.0 / np.pi #c: cm-3 * photons*cm2/keVphot * NTel/cm2/s/keVNT * dkeVNT
      #  print(Q_E.shape)   **** This is equation 2.1 of Holman et al 2011 **** for thick use equation 2.6  ****
    #    e = 0
   #     for e in range(len(E)):
   #         inds = np.all([E < E[e], E > eps_arr[ee]], axis=0)
   #         Q2_E[e] = np.trapz( Q_E[ee, inds] * E[inds], x = E[inds])
    #        j_thick[ee] = j_thick[ee] + np.trapz( Q2_E[e] * F_NN[ee:,], x=E[ee:,]) / 1.1 / (3e-18)  # where K = 3e-18 * LAM_ee / 23 in keV cm2

    ee = 0

    if calc_cttm == 1:
        for ee in range(len(eps_arr)):
            e = 0
            realE = np.all([E > eps_arr[ee]],axis=0)
            Q_E[ee,realE] = 1.4 * 7.90e-25 / eps_arr[ee] / E[realE] * np.log( (1.0 + np.sqrt(1.0 - eps_arr[ee]/E[realE])) / (1.0 - np.sqrt(1.0 - eps_arr[ee]/E[realE])) ) 
            for e in range(len(E)):

                inds = np.all([E < E[e], E > eps_arr[ee]], axis=0)
                Integ1[e] = np.trapz( Q_E[ee, inds] * E[inds], x = E[inds])

            I_thick[ee] = np.trapz( Integ1[ee:,] * F_NN[ee:,], x=E[ee:,]) / 1.1 / (3e-18)  # where K = 3e-18 * LAM_ee / 23 in keV cm2
            # I_thick is in photons/cm2/s/keV

                        

        
    Nphot_HXR = np.trapz(j_eps, x=eps_arr)  # photons/s/cm3/sr
    #c: j_eps is photons/s/cm3/sr/keVphot,  j_lam is ergs/s/cm3/sr/Ang for thin target
    wl_ang = cnst.CC * cnst.HPL / (eps_arr * keV2erg) * 1e8  # angstroms
    j_lam = j_eps * cnst.CC * cnst.HPL / (wl_ang/1e8) * cnst.CC * cnst.HPL / (wl_ang/1e8)**2 / 1e8 * erg2keV
    I_lam_thick = I_thick * cnst.CC * cnst.HPL / (wl_ang/1e8) * cnst.CC * cnst.HPL / (wl_ang/1e8)**2 / 1e8 * erg2keV

    gyrosy_Hz_grid = 10**(np.linspace(0, 3, 100)) * 1e9 # 1 GHz to 1 THz
    frac_mirrored = 1.0
    gyrosy_nub, gyrosy_wl, gyrosy_jlam, Nphot_radio, nu_peak, gyrosy_jHz, gyro_tau, gyro_IHz, Teff_gyro = gyrosynchPL(gyrosy_Hz_grid, Bgauss=Bg, in_theta=45.,delta=delta, N10keV=NN_iso * frac_mirrored, L = radio_L)
    sy_nub, sy_wl, sy_jlam, Nphot_sy, sy_nu_peak, sy_jHz,sy_tau, sy_IHz = synchPL(gyrosy_Hz_grid, Bgauss=Bg, in_theta=45.,delta=delta, N10MeV=np.interp(10000,E_keV,NN_iso_E) * frac_mirrored,Eo=10000.0)

    #c:  in reality, some fraction of the F13 beam precipitates
    #opt_lam_grid = np.load('gen_cont_wave.ccx.pysav.npy')
    #opt_nu_grid = cnst.CCANG/opt_lam_grid
    #Ogyrosy_nub, Ogyrosy_wl, Ogyrosy_jlam, ONphot_radio, Onu_peak, Ogyrosy_jHz,xxtau, xxIHz, xxTeff = gyrosynchPL(opt_nu_grid, Bgauss=Bg, in_theta=45.,delta=delta, N10keV=NN_iso * frac_mirrored, L = radio_L)
    #Osy_nub, Osy_wl, Osy_jlam, ONphot_sy, Osy_nu_peak, Osy_jHz, xytau, xyIHz = synchPL(opt_nu_grid, Bgauss=Bg, in_theta=45.,delta=delta, N10MeV=np.interp(10000,E_keV,NN_iso_E)* frac_mirrored,Eo=10000.0)

    max_E_keV = np.max(E)
    
    #hxr_dict = {'wl_ang':wl_ang, 'eps_arr':eps_arr, 'j_keV':j_eps, 'E_keV':E_keV, 'j_lam':j_lam, 'I_keV_thick':I_thick, 'I_lam_thick':I_lam_thick,'ElNum':F, 'ElDens':NN_iso,'ElDensE':NN_iso_E, 'E10':E10, 'Eldens10':NN_iso10, 'Beam_E':E, 'Beam_Dist':F_NN, 'Nphot_HXR':Nphot_HXR, 'Q':Q_E, 'nuB_Hz':gyrosy_nub, 'gyrosy_jlam':gyrosy_jlam, 'gyrosy_jHz':gyrosy_jHz, 'gyrosy_mm':gyrosy_wl/1e7, 'gyrosy_nu':gyrosy_Hz_grid, 'Nphot_gyrosy':Nphot_radio, 'gyrosy_peakfreq':nu_peak, 'sy_peakfreq':sy_nu_peak, 'sy_jlam':sy_jlam, 'sy_jHz':sy_jHz, 'Nphot_sy':Nphot_sy, 'max_E_keV':max_E_keV, 'opt_lam_Ang':opt_lam_grid, 'OPT_gyrosy_jlam':Ogyrosy_jlam, 'OPT_sy_jlam':Osy_jlam, 'gyrosy_IHz':gyro_IHz,'gyrosy_tau':gyro_tau,'sy_tau':sy_tau, 'sy_IHz':sy_IHz}  # Q is cross section [photenergy, elecenergy] in units of cm2/photkev
    hxr_dict = {'wl_ang':wl_ang, 'eps_arr':eps_arr, 'j_keV':j_eps, 'E_keV':E_keV,
                'j_lam':j_lam, 'I_keV_thick':I_thick, 'I_lam_thick':I_lam_thick,
                'ElNum':F, 'ElDens':NN_iso,'ElDensE':NN_iso_E, 'E10':E10, 
                'Eldens10':NN_iso10, 'Beam_E':E, 'Beam_Dist':F_NN,
                'Nphot_HXR':Nphot_HXR, 'Q':Q_E, 'nuB_Hz':gyrosy_nub, 
                'gyrosy_jlam':gyrosy_jlam, 'gyrosy_jHz':gyrosy_jHz, 
                'gyrosy_mm':gyrosy_wl/1e7, 'gyrosy_nu':gyrosy_Hz_grid, 
                'Nphot_gyrosy':Nphot_radio, 'gyrosy_peakfreq':nu_peak, 
                'sy_peakfreq':sy_nu_peak, 'sy_jlam':sy_jlam, 'sy_jHz':sy_jHz,
                'Nphot_sy':Nphot_sy, 'max_E_keV':max_E_keV,  
                'gyrosy_IHz':gyro_IHz,'gyrosy_tau':gyro_tau,'sy_tau':sy_tau, 
                'sy_IHz':sy_IHz}  # Q is cross section [photenergy, elecenergy] in units of cm2/photkev
    return hxr_dict

def power_law(E_dist, rel_vel,  F_dist_in, n_dist_in, PASIG=15.0, PASIG_METHOD='theta', fine_mu = True):
    keV2erg = 1.602e-9
    # see pasig.ipynb in ccmodel/final.v6/ for some tests and how to create symmetric outputs 
    if PASIG_METHOD == 'theta' or PASIG_METHOD == 'mu' or PASIG_METHOD=='iso':
        if fine_mu:
            nmu_x=60  # number of mu will actually be this number divided by 2
            mu_x = np.zeros(nmu_x)
            for i in range(1,int(nmu_x/2)):
                mu_x[i-1] = np.cos((i)/(nmu_x/2.0) * np.pi/2)**.85

            good_mu = (mu_x > 0.0)
            mu_xx = mu_x[good_mu]
            theta_grid = np.arccos(mu_xx) * 180.0 / np.pi

        else:
            theta_grid_a = np.arange(0.1,10.0,0.1)
        
            theta_grid_b = np.arange(10.0,90.0,1.0) # 0-90 degrees.

            theta_grid = np.append(theta_grid_a, theta_grid_b)
        ###
       # theta_grid = np.arange(1.0,90.0,1.0)
        
        theta_grid_rad = theta_grid /180.0 * np.pi
        mu_grid = np.cos(theta_grid_rad)
        mu_grid1 = prep_pmesh(mu_grid)
        dmu_grid = np.abs(mu_grid1[1:-1] - mu_grid1[0:-2])
        dmu_grid = np.append(dmu_grid, dmu_grid[-1])

        E_dist1 = prep_pmesh(E_dist)
        dE_grid = E_dist1[1:-1] - E_dist1[0:-2]
        dE_grid = np.append(dE_grid, dE_grid[-1])

        n_energy_bins = len(E_dist)
        n_theta_bins = len(theta_grid_rad)
        F_dist_out = np.zeros((n_energy_bins, n_theta_bins))

        #F_dist_in should be in erg/cm2/s/keV, E_dist in keV
        F_nrm = np.trapz(F_dist_in, x=E_dist)  # erg/s/cm2
        print('Input erg/s/cm2 = ',F_nrm)

        if PASIG_METHOD == 'theta':
            theta0 = PASIG
            theta0_rad = theta0 / 180.0 * np.pi
            M =  np.exp( -1.0 * (theta_grid_rad)**2 / (2.0 * theta0_rad**2) ) / np.sqrt(2.0 * np.pi) / theta0_rad  # Mtheta
        if PASIG_METHOD == 'mu':
            mu0 = PASIG / np.sqrt(2.0)
            M =  np.exp( -1.0 * (mu_grid - 1.0)**2 / (2.0 * mu0**2) ) / np.sqrt(2.0 * np.pi) / mu0  # M mu from Allred et al. 2015 (fp_solver.f)
        if PASIG_METHOD == 'iso':
            M = mu_grid/mu_grid
        

        rel_vel_2d = np.broadcast_to(rel_vel, (n_theta_bins, n_energy_bins))
        rel_vel_2D = np.transpose(rel_vel_2d)
        n_dist_in_2d = np.broadcast_to(n_dist_in, (n_theta_bins,n_energy_bins))
        n_dist_in_2D = np.transpose(n_dist_in_2d)
        Itemp = np.broadcast_to(M, (n_energy_bins, n_theta_bins)) * rel_vel_2D * n_dist_in_2D #  number/cm2/s/keV as a funtion of theta in radians.
        nrm = 0.0
        nrm2 = 0.0

     #   print(len(E_dist), len(mu_grid), Itemp.shape, dE_grid.shape, dmu_grid.shape, rel_vel.shape, mu_grid.shape, E_dist.shape, theta_grid_rad.shape)
        for nEE in range(len(E_dist)):
            for nUU in range(len(mu_grid)):                
                nrm += 2.0 * np.pi * Itemp[nEE, nUU] * dE_grid[nEE] * dmu_grid[nUU] * mu_grid[nUU] * E_dist[nEE] * keV2erg
                nrm2 += 2.0 * np.pi * Itemp[nEE, nUU] * dE_grid[nEE] * dmu_grid[nUU] *  mu_grid[nUU]


        I = Itemp * F_nrm / nrm # properly normalized intensity:  number/s/cm2/steradian/keV, such that when integrated over forward hemisphere (and then over energy) will give desired erg/s/cm2 such as F13 or whatever..   integral( E * I * mu * domega) = flux and integral(flux * dE) = F#.


        E_grid_2d = np.broadcast_to(E_dist, (n_theta_bins, n_energy_bins))
        E_grid_2D = np.transpose(E_grid_2d)
        Intensity = I * (E_grid_2D * keV2erg) # erg/s/cm2/sr/keV
        
        check1= 0.0
        for nEE in range(len(E_dist)):
            for nUU in range(len(mu_grid)):   
                check1 += 2.0 * np.pi * I[nEE, nUU] * dE_grid[nEE] * dmu_grid[nUU]  * mu_grid[nUU]

        check2= 0.0
        for nEE in range(len(E_dist)):
            for nUU in range(len(mu_grid)):   
                check2 += 2.0 * np.pi * Intensity[nEE, nUU] * dE_grid[nEE] * dmu_grid[nUU] *  mu_grid[nUU] 

        check3= 0.0
        for nEE in range(len(E_dist)):
            for nUU in range(len(mu_grid)):   
                check3 += 2.0 * np.pi * I[nEE, nUU] * dE_grid[nEE] * dmu_grid[nUU] / rel_vel[nEE]  


    print(check1, check2, check3)
    return E_dist, dE_grid, F_dist_in, mu_grid, dmu_grid, theta_grid, I, Intensity, check1, check2, check3 #  keV, erg/s/cm2/keV (nE), mu_grid (nMu), number/s/cm2/sr/keV, erg/s/cm2/sr/keV, number/s/cm2, erg/s/cm2, number/cm3
        
        
      #  nrm_over_E = np.trapz(2.0 * np.pi * Itemp * mu_grid * rel_vel
        
      #  nrm_gauss = np.trapz(gaubeam * sin(theta_grid_rad), x = theta_grid_rad)  # integrate from -90 to 90.

      #  padist = gaubeam / (nrm  * 2.0 * np.pi)
# check that spherical coordinate integral over padist is 1.0.

    # integral of Mtheta * domega = 1.0 b/c Mtheta must be in units of per steradian.
   #     Mtheta = gaubeam
    # now solve for N0 given F0 and then plug that into Joel's equation for F0 to check for consistency.
    # for a given N0, solve for energy flux

    
  ####  F0 = N0 * (delta-1.0) / (Ec * keV2erg) * Mtheta * (E_grid / Ec)**(-1.0 * delta)  # integral over all energies is N0 * Mtheta.
    # from joel:  F0 in paper is particles / cm^2 /s^1 / keV / sr
    # M0 sintheta dtheta dphi = 4pi steradians, this can't be right.

    # then take the integral with cos theta to get the flux (make sure to integrate over both hemis for the flux because of the 4pi in equation for A!)
    #  check fp to see if this is true because it's only gaussian in forward hemisphere.

    # can do the same thing with emergent radiative intensity.
