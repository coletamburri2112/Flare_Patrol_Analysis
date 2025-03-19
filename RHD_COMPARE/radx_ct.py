import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#import astropy
#from scipy import interpolate
import scipy.io
import radxglobals as cnst
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

def format_coord(x, y):
    global X, Y, ZVALS
    xarr=X
    yarr=Y
    colx = findind(xarr,x)
    rowy = findind(yarr,y)
    zval = ZVALS[rowy, colx]
    return 'x=%1.4f, y=%1.4f, indx = %1i, indy=%4i, val=%1.4e' % (x, y, colx, rowy, zval)

def format_coordxy(x, y):
    global XX, YY, ZZ
    xarr=XX
    yarr=YY
    zarr=ZZ
    colx = findind(xarr,x)
    #rowy = findind(yarr,y)
    #rowy2 = findind(zarr,
    zval = zarr[colx]
    yval = yarr[colx]
    xval = xarr[colx]
    return 'x=%1.4f, y1=%1.2e, y2=%1.2e, indx = %1i' % (xval, zval, yval, colx)


def prep_pmesh(z):
    # z is an irregular grid and must be at cell boundaries for pcolormesh (therefore make an array that is ndep + 1 dimensions.)
    ndep = len(z)
    midz = (z[1:len(z)] + z[0:len(z)-1])/2.
    newz = np.insert(midz, 0, z[0] + (z[0]-midz[0]))
    ndep2=len(newz)
    z_bdry = np.append(newz, z[ndep-1] + (z[ndep-1]-midz[ndep-2]))
    return z_bdry

def findind(array,value):
    ''' closest_index = findind(array, value_of_interest);  finds index of array where array is closest to the given value.'''
    idx = (np.abs(array-value)).argmin()
    return idx

def fi(array,value):
    '''closest_index = fi(array, value_of_interest); finds index of array where array is closest to the given value.'''
    idx = (np.abs(array-value)).argmin()
    return idx

def main():
    return None

def help():
    names = ['Height','Col Mass', 'Temp', 'Density', 'Nel', 'Beam Heating', 'Pressure', 'Velocity']
    info = ['z1t/1e5 [km]', 'log10 cmass1t [g/cm2]', 'tg1t [K]', 'd1t [g/cm3]', 'ne1t [electrons/cm3]', 'bheat1t [erg/s/cm3]', 'pg1t [dyn/cm2]','vz1t/1e5 [km/s], positive=up']
    print('The options for rad2plot are the following: ')
    for nn in range(len(names)):
        print(names[nn],':', info[nn])
    return None
def rad2plot(atmos,x1in, y1in, y2in, time = 0.0, xlim=[-100,1500], y1lim=[3,6], y2lim=[3,6], user_figsize=(8,6), user_fontsize=16, oplot_t0=False,psym=False, y1log = False, y2log = False,savefig=False,plot_filename='rad2plot',y2_other=0):
    if x1in == 'Height':
        x1 = atmos.z1t/1e5
        x1in = x1in + ' (km)'
    if x1in == 'Col Mass':
        x1 = np.log10(atmos.cmass1t)
        x1in = x1in + r' (g cm$^{-2}$)'
        x1in = r'log$_{10}$ '+x1in
        
    if y1in == 'Temp':
        y1 = atmos.tg1t
        y1in = y1in +' (K)'
    if y1in == 'Density':
        y1 = atmos.d1t
        y1in = y1in+r' (g cm$^{-3}$)'
    if y1in == 'Nel':
        y1 = atmos.ne1t
        y1in = r'N$_e$ (cm$^{-3}$)'
    if y1in == 'Beam Heating':
        y1 = atmos.bheat1t
        y1in = y1in + r' (erg s$^{-1}$ cm$^{-3}$)'
    if y1in == 'Pressure':
        y1 = atmos.pg1t
        y1in = y1in + r' (dyn cm$^{-2}$)'
    if y1in == 'Velocity':
        y1 = atmos.vz1t/1e5
        y1in = y1in+r' (km s$^{-1}$)'

    if y2in == 'Temp':
        y2 = atmos.tg1t
        y2in = y2in +' (K)'
    if y2in == 'Density':
        y2 = atmos.d1t
        y2in = y2in+r' (g cm$^{-3}$)'
    if y2in == 'Nel':
        y2 = atmos.ne1t
        y2in = r'N$_e$ (cm$^{-3}$)'
    if y2in == 'Beam Heating':
        y2 = atmos.bheat1t
        y2in = y2in + r' (erg s$^{-1}$ cm$^{-3}$)'
    if y2in == 'Pressure':
        y2 = atmos.pg1t
        y2in = y2in + r' (dyn cm$^{-2}$)'
    if y2in == 'Velocity':
        y2 = atmos.vz1t/1e5
        y2in = y2in+r' (km s$^{-1}$)'
    if y2in == 'user':
        y2 = y2_other
        
    if y1log:
        y1 = np.log10(y1)
        y1in = r'log$_{10}$ '+y1in
    if y2log:
        y2 = np.log10(y2)
        y2in= r'log$_{10}$ '+y2in
    
    timet = atmos.timet
    brightcol = color_bright()
    rnbw_col = color_rainbow14()
    plt.rc('font',size=user_fontsize)
    f, ax1 = plt.subplots(figsize=user_figsize)
    indt1 = findind(timet, time)
    utime = timet[indt1]
    print('time = ',timet[indt1])
    print('min max of y1:',np.min(y1[:,indt1]), np.max(y1[:,indt1]))
    print('min max of y2:',np.min(y2[:,indt1]), np.max(y2[:,indt1]))

    if psym:
        ax1.plot(x1[:,indt1], y1[:,indt1], ls='solid',color='k',lw=2,marker='+')
    else:
        ax1.plot(x1[:,indt1], y1[:,indt1], ls='solid',color='k',lw=2)
    if oplot_t0 == True:
        ax1.plot(x1[:,0], y1[:,0], ls='dotted',color='k')
        print('Dotted linestyles are the t=0s atmospheric parameters.')
    ax1.set_ylim(y1lim)
    ax1.set_xlim(xlim)
    ax1.set_xlabel(x1in)
    ax1.set_ylabel(y1in)
    ax2 = ax1.twinx()
    if psym:
        ax2.plot(x1[:,indt1], y2[:,indt1],ls='dashed',color=rnbw_col[13],lw=2,marker='+')
    else:
        ax2.plot(x1[:,indt1], y2[:,indt1],ls='dashed',color=rnbw_col[13],lw=2)
    if oplot_t0 == True:
        ax2.plot(x1[:,0], y2[:,0],ls='dashdot',color=rnbw_col[13])
    ax2.set_ylabel(y2in,color=rnbw_col[13])
    ax2.set_ylim(y2lim)
    
    ax2.tick_params(axis='y', colors=rnbw_col[13])
    ax2.spines['right'].set_color(rnbw_col[13])
    ax2.set_title('t = '+str('{0:.2f}'.format(utime))+' s')
    global XX
    global YY
    global ZZ
    XX = x1[:,indt1]
    YY = y2[:,indt1]
    ZZ = y1[:,indt1]
    ax2.format_coord = format_coordxy
    if savefig:
        plt.savefig(plot_filename+'.pdf')
        print('created plot: ',plot_filename+'.pdf')
    return None

class modelclass:
    pass

def load_atmos(H_2 = False):
    model_dir = np.load('model_dir.tmp.npy')
    model_file = np.load('model_file.tmp.npy')
    radynvar = scipy.io.readsav(str(model_dir)+str(model_file))
    atmos = modelclass()
    atmos.z1t = radynvar.z1t.transpose()  # have to transpose all this stuff because scipy.io.readsav transposes.   atmos.z1t[ndep, ntime]
    atmos.timet = radynvar.timet
    atmos.vz1t = radynvar.vz1t.transpose()
    atmos.d1t = radynvar.d1t.transpose()
    atmos.zmu = radynvar.zmu
    atmos.dzt = radynvar.dzt.transpose()
    atmos.pg1t = radynvar.pg1t.transpose()
    atmos.ne1t = radynvar.NE1T.transpose()
    atmos.taut = radynvar.taut.transpose()
    atmos.z1t = radynvar.Z1T.transpose()
    atmos.tg1t = radynvar.TG1T.transpose()
    atmos.n1t = radynvar.N1T.transpose()
    atmos.totnt = radynvar.TOTNT.transpose()
    try:
        atmos.nstart = radynvar.NSTART.transpose()
        atmos.c1t = radynvar.CT.transpose()  # collisional rates.
        atmos.rijt = radynvar.RIJT.transpose()
        atmos.rjit = radynvar.RJIT.transpose()
        atmos.f20t = radynvar.F20T.transpose()

    except:
        print('Could not read in nstart, rijt, c1t.')
    atmos.cmass1t = radynvar.CMASS1T.transpose() 
    atmos.bheat1t = radynvar.bheat1t.transpose()
    atmos.cont = radynvar.cont
    atmos.irad = radynvar.irad
    atmos.jrad = radynvar.jrad
    atmos.alamb = radynvar.alamb
    atmos.ielrad = radynvar.ielrad
    atmos.atomid = radynvar.atomid
    atmos.label = radynvar.label
    atmos.ion = radynvar.ion.transpose()
    atmos.g = radynvar.g.transpose()
    atmos.grph = radynvar.GRPH
    if H_2:
        atmos.H_2 = (atmos.d1t / atmos.grph - atmos.totnt[:,0,:])*0.5 # if H_2 populations were included in LTE chemical equilibrium in radyn (m dwarf only)
        print('Read in H_2 population densities.')
    #isel = np.all([irad == 2],axis=0)
    #jsel = np.all([jrad == 2],axis=0)
    return atmos


def printatmos(atmos, time = 0.0, height_index = -99):
    timet = atmos.timet
    indt1 = findind(timet, time)
    if height_index < 0:  # tau500 = 1 closest index
        height_index = findind(atmos.taut[:,indt1], 1.0)
        print('Using height index = ',height_index, ' for tau500=',atmos.taut[height_index,indt1])
    print('******')
    print('Params at time [s] = ',timet[indt1])
    print('z [km] = {0:.3f}'.format(atmos.z1t[height_index, indt1]/1e5))
    print('log10 m [g cm-2] = {0:.3e}'.format( np.log10(atmos.cmass1t[height_index, indt1])))
    tg = atmos.tg1t[height_index, indt1]
    nelg = atmos.ne1t[height_index, indt1]
    dg = atmos.d1t[height_index, indt1]
    heIIIfrac = atmos.n1t[height_index, 8, 2, indt1] / atmos.totnt[height_index,2,indt1]
    print('Temp [K] = {0:.3f}, N_electron [cm-3] = {1:.3e}, rho [g cm-3] = {2:0.3e}'.format(tg, nelg, dg))
    print('He III fraction = {0:.6e}'.format(heIIIfrac))
    print('******')

    return None


def interpatmos(atmos, time = 0.0, height=0.0, tau = False):  # height in km if tau is false
    timet = atmos.timet
    indt1 = findind(timet, time)
    if tau:
        ival = 1.0
        iarr = atmos.taut[:,indt1]
        t_interp = np.interp(ival,iarr,atmos.tg1t[:,indt1])
        z_interp = np.interp(ival,iarr, atmos.z1t[:,indt1]/1e5)
        m_interp = np.interp(ival,iarr, np.log10(atmos.cmass1t[:,indt1]))
        ne_interp = np.interp(ival,iarr,atmos.ne1t[:,indt1])
        d_interp = np.interp(ival,iarr,atmos.d1t[:,indt1])
    if not tau:
        ival = height
        iarr = np.flip(atmos.z1t[:,indt1]/1e5)
        t_interp = np.interp(ival,iarr,np.flip(atmos.tg1t[:,indt1]))
        z_interp = np.interp(ival,iarr, np.flip(atmos.z1t[:,indt1]/1e5))
        m_interp = np.interp(ival,iarr, np.flip(np.log10(atmos.cmass1t[:,indt1])))
        ne_interp = np.interp(ival,iarr,np.flip(atmos.ne1t[:,indt1]))
        d_interp = np.interp(ival,iarr,np.flip(atmos.d1t[:,indt1]))
    print('*******')
    print('Params at tau500 = 1:')
    print('--------------------')
    print('time [s] = ',timet[indt1])
    print('z [km] = {0:.3f}'.format(z_interp))
    print('log10 m [g cm-2] = {0:.3e}'.format(m_interp ))
    print('Temp [K] = {0:.3f}, N_electron [cm-3] = {1:.3e}, rho [g cm-3] = {2:0.3e}'.format(t_interp, ne_interp, d_interp))
    print('******')
    return None




def lineguide(modx):
    cont = modx.cont
    irad = modx.irad
    jrad = modx.jrad
    alamb = modx.alamb
    ielrad = modx.ielrad
    atomid = modx.atomid
    ion = modx.ion
    gi = modx.g
    print('type    i     j  lam_cen (vac)  EL+Stage     kr(ignore)')
    for i in range(len(cont)):
        kr = i
        iel = ielrad[i] - 1
        if ielrad[i] == 1:
            ielradid = 'H'
        if ielrad[i] == 2:
            ielradid = 'Ca II'
        if ielrad[i] == 3:
            ielradid = 'He'
        if cont[i] != 0:
            zindype = 'b-f'
        if cont[i] == 0:
            zindype = 'b-b'
        print('{0:s} {1:5.0f} {2:5.0f} {3:<15.2f} {4:2s} {5:<5.0f} {6:10.0f}'.format(zindype, irad[i], jrad[i], alamb[i], ielradid, ion[irad[i]-1,ielrad[i]-1], kr))
    
    return None


def vactoair(wave_vac):

    sigma2 = (1e4/(wave_vac) )**2   
    fact = 1.0 +  5.792105E-2/(238.0185E0 - sigma2) + 1.67917E-3/( 57.362E0 - sigma2)
    wave_air = wave_vac/fact

    return wave_air



# Good color schemes to use from https://www.sron.nl/~pault/
def color_map(umap = 'rnbw'):
#    ''' user_cmap = mf.color_map(umap='rnbw') where umap can be burd, burd_flip, or bryl'''
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib import cm
    if umap == 'rnbw':  # this is rainbow34 aka rainbow_WhBr from Figure 20 of Paul Tol's website for interpolating.
        print('Brown to White rainbow.')
        clrs = ['#E8ECFB', '#DDD8EF', '#D1C1E1', '#C3A8D1', '#B58FC2','#A778B4','#9B62A7', '#8C4E99', '#6F4C9B', '#6059A9',  '#5568B8', '#4E79C5', '#4D8AC6', '#4E96BC', '#549EB3', '#59A5A9', '#60AB9E', '#69B190', '#77B77D', '#8CBC68',  '#A6BE54', '#BEBC48', '#D1B541', '#DDAA3C', '#E49C39', '#E78C35', '#E67932', '#E4632D', '#DF4828', '#DA2222', '#B8221E', '#95211B', '#721E17', '#521A13']
        cmap_name = 'rainbow_brwh'
        usermap = LinearSegmentedColormap.from_list(cmap_name, np.flip(clrs), N=500)
      #  usermap.set_bad('#666666')
        usermap.set_bad('#521A13')
    if umap == 'rnbw_flip':  # this is rainbow34 aka rainbow_WhBr from Figure 20 of Paul Tol's website for interpolating.
        print('Brown to White rainbow.')
        clrs = ['#E8ECFB', '#DDD8EF', '#D1C1E1', '#C3A8D1', '#B58FC2','#A778B4','#9B62A7', '#8C4E99', '#6F4C9B', '#6059A9',  '#5568B8', '#4E79C5', '#4D8AC6', '#4E96BC', '#549EB3', '#59A5A9', '#60AB9E', '#69B190', '#77B77D', '#8CBC68',  '#A6BE54', '#BEBC48', '#D1B541', '#DDAA3C', '#E49C39', '#E78C35', '#E67932', '#E4632D', '#DF4828', '#DA2222', '#B8221E', '#95211B', '#721E17', '#521A13']
        cmap_name = 'rainbow_brwh'
        usermap = LinearSegmentedColormap.from_list(cmap_name, clrs, N=500)
    elif umap == 'burd':
        BuRd = ["#2166AC", "#4393C3", "#92C5DE", "#D1E5F0","#F7F7F7", "#FDDBC7","#F4A582", "#D6604D", "#B2182B"]  # bad = 255,238,153 = FFEE99
        cmap_name = 'BuRd' 
        usermap = LinearSegmentedColormap.from_list(cmap_name, np.flip(BuRd), N=100)
    elif umap == 'burd_flip':
        BuRd = ["#2166AC", "#4393C3", "#92C5DE", "#D1E5F0","#F7F7F7", "#FDDBC7","#F4A582", "#D6604D", "#B2182B"]  # bad = 255,238,153 = FFEE99
        cmap_name = 'BuRd_Flipped' 
        usermap = LinearSegmentedColormap.from_list(cmap_name, BuRd, N=100)
    elif umap == 'bryl':
        clrs_ylbr = ['#FFFFE5', '#FFF7BC','#FEE391','#FEC44F','#FB9A29','#EC7014','#CC4C02', '#993404','#662506']
        cmap_name = 'ylbr'
        usermap = LinearSegmentedColormap.from_list(cmap_name, np.flip(clrs_ylbr), N=500)
    else:
        print( ' umap can be rnbw, burd, burd_flip, or bryl')
        
    return usermap
def color_rainbow14(printc = 'no'):
    ''' This is rainbow14 plus grey as last entry, Figure 18 top panel of Paul Tol's website.  color_rainbow(printc = no or yes)'''
    rainbow = [(209,187,215), (174,118,163), (136,46,114), (25,101,176), (82,137,199), (123,175,222), (77,178,101), (144,201,135), (202, 224, 171), (247, 240, 86), (246,193, 65), (241,147,45), (232, 96,28), (220, 5,12), (119, 119, 119)]
    labels=['ltpurple0', 'medpurple1','darkpurple2', 'darkblue3','medblue4', 'lightblue5', 'darkgreen6','medgreen7', 'ltgreen8','yellow9','ltorange10','medorange11', 'dkorange12', 'red13', 'grey14']
    for i in range(len(rainbow)):    
        r, g, b = rainbow[i]    
        rainbow[i] = (r / 255., g / 255., b / 255.)
        if printc == 'yes' or printc =='y':
            print(i, labels[i])
    return rainbow


def color_bright(printc='no'):
    ''' color_bright(printc = no or yes) '''
    bright = [(68,119,170), (102,204,238), (34, 136, 51), (204,187,68), (238,102,119), (170,51,119), (187,187,187)]   
    labels=['blue' ,'cyan', 'green', 'yellow','red','purple', 'light grey']
    for i in range(len(bright)):    
        r, g, b = bright[i]    
        bright[i] = (r / 255., g / 255., b / 255.)
        if printc == 'yes' or printc =='y':
            print(i, labels[i])
    return bright

   
        
def ci_image1(atmos,line='Ha',mu=0.95,time=0.,xlim=[6560,6564],ylim=[-100,1500],
              ci_log=False,vmax=0.5,vmin=0.5*1e-5,user_cmap='gray', 
              oplot_Ilam = True, Ilam_color='white', oplot_vel=True,
              vel_color='gray',oplot_legend2=True,oplot_legend1=True,
              savefig=False,plot_filename='contrib_fn_plot',user_figsize=(6,6),
              user_fontsize=14,ci_prime=[-99,-99],ci_zlim=[-200,-200],
              ciprime_color='white',oplot_tau1=False, tau_color='k',
              tau_2d = False,tlim0=0,tlim1=30,src_2d=False,source_ha=False,
              oplot_t0=True, oplot_annotate=False,flip_wl=True,vel_ls='dashdot'):
    '''  Example:
    your_python_dictionary = radx.ci_image1(atmos,time=18,vmin=-2,vmax=1,xlim=[6560,6566],\
        user_cmap=your_color_map,savefig=True,user_figsize=(8,6),ci_log=True,user_fontsize=14,\
        vel_color='#BBBBBB',oplot_tau1=1) '''
    global lam_rest
    model_dir = np.load('model_dir.tmp.npy')
    model_contrib = np.load('model_contrib.tmp.npy')
    contribf = scipy.io.readsav(str(model_dir)+str(model_contrib),verbose=False,python_dict=True)
    model_file = np.load('model_file.tmp.npy')
    cf = contribf['lcontribf'][0]
  #  cf.dtype
    col_rnbw = color_rainbow14()
    if mu < 0.51:
        print('Are you sure you want to look at a mu-value <= 0.5 for a plane parallel flare atmosphere?')
    muind = findind(atmos.zmu, mu)
    user_mu = atmos.zmu[muind]
    if line == 'Ha' and mu==0.95: 
        cfline = cf.LHAINTT95
        plt_label = r'H$\alpha$'
        lam_rest = atmos.alamb[2]
    if line == 'Ha' and mu==0.77: 
        cfline = cf.LHAINTT77
        plt_label = r'H$\alpha$'
        lam_rest = atmos.alamb[2]
        
    if line == 'HeII' and mu==0.95:
        cfline = cf.LHEII304INTT95
        plt_label = r'H$\alpha$'
        lam_rest = atmos.alamb[31]
    if line == 'HeII' and mu==0.77:
        cfline = cf.LHEII304INTT77
        plt_label = r'H$\alpha$'
        lam_rest = atmos.alamb[31]
        
    if line == 'Hg' and mu ==0.95:
        cfline = cf.LHGINTT95
        plt_label=r'H$\gamma$'
        lam_rest = atmos.alamb[7]
        # need to finish
    if line == 'Hg' and mu ==0.77:
        cfline = cf.LHGINTT77
        plt_label=r'H$\gamma$'
        lam_rest = atmos.alamb[7]
        # need to finish

    if line == 'Ca II K' and mu ==0.95:
        cfline = cf.LCAIIKINTT95
        plt_label = 'Ca II K'
        lam_rest = atmos.alamb[17]
    if line == 'Ca II 8542' and mu ==0.95:
        cfline = cf.LCAII8542INTT95
        plt_label = 'Ca II 8542'
        lam_rest = atmos.alamb[20]
    if source_ha == True:
        cfline_ha = cf.LHAINTT95
        # need to finish

    if lam_rest > 3000.:
        lam_rest = vactoair(lam_rest)
        
    #print('Plotting contribution function (erg/s/cm3/sr/Ang) for ',plt_label,' at mu= {:0.3f}'.format(user_mu), 'at time = {:0.3f} '.format(time), 'sec for model',str(model_file),'.')
    # to do:  put in other lines that adam has stored.
    #print(cfha.dtype)   this prints the available arrays 
    contribline = cfline['contribf']
    contriblinewl = cfline['lam']
    c_I_line = contribline[0].transpose()  # at mu = 0.95
    c_I_line_wl = contriblinewl[0].transpose()


   # contribline_ha = cfline_ha['contribf']
   # contriblinewl_ha = cfline_ha['lam']
    
    tind = findind(atmos.timet, time)
    if ylim[0] < 0:
        ylim[0] = np.min(atmos.z1t[:,tind]/1e5)
    zprep = prep_pmesh(atmos.z1t[:,tind]/1e5)
    zprep2 = atmos.z1t[:,tind]/1e5
    wlprep2 = c_I_line_wl[:,tind]
    wlprep = prep_pmesh(c_I_line_wl[:,tind])
    wl2d, z2d = np.meshgrid(wlprep,zprep)
    vz2d, z2dx = np.meshgrid((lam_rest-wlprep)/lam_rest*2.998e5,zprep)
    vzprep2 = (lam_rest-wlprep2)/lam_rest*2.998e5
    user_time = atmos.timet[tind]
    plt.rc('font',size=user_fontsize)
    f, (ax1,ax4) = plt.subplots(ncols=2,figsize=user_figsize, sharey='row',
                               gridspec_kw={'width_ratios': [4,1]})
#fig.set_tight_layout({'rect': [0, 0, 1, 0.95], 'pad': 0, 'h_pad': 0})
    if user_cmap != 'gray' and user_cmap != 'bone':
        xmap = color_map(umap=user_cmap)
    else:
        xmap = user_cmap

    if ci_log == 0:
        ax1.pcolormesh(wl2d, z2d, c_I_line[:,:,tind],vmax=vmax,vmin=vmin,cmap = xmap,rasterized=True)
    else:
        ci0 = np.log10(c_I_line[:,:,tind].squeeze())
        if tau_2d == True:
            tauall = cfline['tau']
            tauline = tauall[0].transpose()
         #   print(tauline.shape, wl2d.shape, z2d.shape)
            if flip_wl:
                ax1.pcolormesh(wl2d, z2d,np.log10(np.flip(tauline[:,:,tind]/mu,axis=1)),vmax=vmax,vmin=vmin,cmap = xmap,rasterized=True)
            else:
                ax1.pcolormesh(wl2d, z2d,np.log10(tauline[:,:,tind]),vmax=vmax,vmin=vmin,cmap = xmap,rasterized=True)
        elif src_2d == True:
            srcall = cfline['src']
            srcline = srcall[0].transpose()
            #print(srcline.shape, wl2d.shape, z2d.shape)
            plvec = np.vectorize(planckfni)
            b1d = np.zeros((len(atmos.tg1t[:,tind]),1))
            b1d[:,0] = plvec(np.median(wlprep), atmos.tg1t[:,tind])
            if flip_wl:
                src_flipped = np.flip(srcline[:,:,tind],axis=1)
            else:
                src_flipped = srcline[:,:,tind]
            dnu2dlam = np.zeros_like(src_flipped) + cnst.CCANG / np.median(wl2d)**2
            src_flipped_dlam = src_flipped * dnu2dlam
            b2d = np.broadcast_to(b1d, src_flipped.shape)
            #print(b1d)
            ax1.pcolormesh(wl2d, z2d,( src_flipped_dlam/ b2d),vmax=vmax,vmin=vmin,cmap =xmap,rasterized=True)
            src_ratio =  src_flipped_dlam/ b2d
        else:
            ax1.pcolormesh(wl2d, z2d,ci0,vmax=vmax,vmin=vmin,cmap = xmap,rasterized=True)
    ax1.set_ylim(ylim)
    ax1.set_xlim(xlim)
    ax1.set_xlabel(r'Wavelength ($\rm{\AA}$)')
    ax1.set_ylabel('Distance [km] from Photosphere')
    if oplot_Ilam:
        contriblineint = cfline['int']
        emerg_intline = contriblineint[0].transpose()
        ax1.plot(c_I_line_wl[:,tind], emerg_intline[:,tind] * ylim[1] / np.max(emerg_intline[:,tind]) * 0.75 ,color=Ilam_color,label=r'I$_{\lambda}$',lw=2.0,marker='+',ms=6.)
        if oplot_legend1:
            ax1.legend(loc='upper left',frameon=False,labelcolor=Ilam_color)
        #print('The maximum emergent I_lam (erg/s/cm2/sr/Ang) for this profile is {:0.3e}'.format(np.max(emerg_intline[:,tind])))

    if oplot_tau1:
        tauall = cfline['tau']
        tauline = tauall[0].transpose()
        shptau = tauline.shape
        tau2deq1 = np.zeros((shptau[1]))
        for tt in range(shptau[1]):
            tau2deq1[tt] = np.interp(1.0, tauline[:,tt,tind] / mu, atmos.z1t[:,tind]/1e5)
        if flip_wl:
            ax1.plot(c_I_line_wl[:,tind], np.flip(tau2deq1),color=tau_color,ls='dashed',label=r'$\tau=1$',lw=2,zorder=22)  # need to flip tau because contribf is flipped in contribfunc.pro
        else:
            ax1.plot(c_I_line_wl[:,tind], tau2deq1,color=tau_color,ls='dashed',label=r'$\tau=1$',lw=2,zorder=22)
        ax1.legend(loc='upper left',frameon=False,labelcolor='white',ncol=2)

        
    if ci_prime[0] > 0 and ci_prime[1] > 0:
        nlam = len(c_I_line_wl[:,tind])
        ciinds = np.all([atmos.z1t[:,tind]/1e5 > ci_zlim[0], atmos.tg1t[:,tind] < 1e5],axis=0)
        nzci = np.count_nonzero(ciinds)
        CIPRIME = np.zeros((nlam,nzci) )
        CIPRIMEZ0 = np.zeros(nlam)
        for ww in range(nlam):
            dzt_sel = atmos.dzt[ciinds,tind]
            c_I_line_sel = c_I_line[ciinds, :, tind]
            z1t_sel = atmos.z1t[ciinds,tind]
            ZPRIME, CIPRIME[ww,:] = np.abs(akcdf_dz(dzt_sel, c_I_line_sel[:,ww],norm=True))
            CIPRIMEZ0[ww] = np.interp(ci_prime[0], CIPRIME[ww,:], z1t_sel)
        ax1.plot(c_I_line_wl[:,tind], CIPRIMEZ0/1e5, ls=(0,(5,1)),color=ciprime_color,lw=0.7)
        
        nlam = len(c_I_line_wl[:,tind])
        ciinds = np.all([atmos.z1t[:,tind]/1e5 > ci_zlim[1], atmos.tg1t[:,tind] < 1e5],axis=0)
        nzci = np.count_nonzero(ciinds)
        CIPRIME = np.zeros((nlam,nzci) )
        CIPRIMEZ1 = np.zeros(nlam)
        for ww in range(nlam):
            dzt_sel = atmos.dzt[ciinds,tind]
            c_I_line_sel = c_I_line[ciinds, :, tind]
            z1t_sel = atmos.z1t[ciinds,tind]
            ZPRIME, CIPRIME[ww,:] = np.abs(akcdf_dz(dzt_sel, c_I_line_sel[:,ww],norm=True))
            CIPRIMEZ1[ww] = np.interp(ci_prime[1], CIPRIME[ww,:], z1t_sel)
        ax1.plot(c_I_line_wl[:,tind], CIPRIMEZ1/1e5, ls=(0,(5,1)),color=ciprime_color,lw=0.7)
    global X 
    X= vzprep2
    global Y
    Y=zprep2
    global ZVALS

    if tau_2d == True:
        if flip_wl:
            ZVALS = np.flip(tauline[:,:,tind] / mu,axis=1)
        else:
            ZVALS = tauline[:,:,tind] / mu
    elif src_2d == True:
        ZVALS = src_ratio
    else:
        ZVALS = c_I_line[:,:,tind]
        
    ax2 = ax1.twiny()
    if oplot_vel:
        ax2.plot(atmos.vz1t[:,tind]/1e5 * mu, atmos.z1t[:,tind]/1e5,ls=vel_ls,color=vel_color,label=r'$v$',lw=1.75)  # ls=(0,(3,1,1,1,1,1))

    dlammin = lam_rest - xlim[0]
    dlamplus = lam_rest - xlim[1]

    ax2.set_xlim(dlammin/lam_rest *2.998e5, dlamplus/lam_rest * 2.998e5)
    ax2.set_xlabel('L.o.S. Gas Velocity (km s$^{-1}$); negative = downward')

    if oplot_annotate:
        ax2.text(-325,100,'Backwarmed Upper',ha='center',va='center',fontsize=10)
        ax2.text(-325,50,'Photosphere',ha='center',va='center',fontsize=10)

        ax2.text(-325,830,'Stationary Chrom',ha='center',va='center',fontsize=10)
        ax2.text(-325,780,'Flare Layers',ha='center',va='center',fontsize=10)
        ax2.text(-325,890,'CC',ha='center',va='center',fontsize=10)
        ax2.text(-325,1000,'Evaporation',ha='center',va='center',fontsize=10)

    if oplot_legend2:
        ax2.legend(loc='upper right',frameon=False,labelcolor=vel_color)
    ax2.format_coord = format_coord
    ax4.plot(atmos.tg1t[:,tind]/1000.0, atmos.z1t[:,tind]/1e5,color='k')
    ax4.set_xlim(tlim0,tlim1)
    bright = color_bright()
    tinc = np.arange(0, tlim1, 5)
    for ti in range(len(tinc)):
        ax4.plot([tinc[ti],tinc[ti]],[-100,1500],color=bright[6],lw=0.7)

    ax4.plot(atmos.tg1t[:,tind]/1000.0, atmos.z1t[:,tind]/1e5,color='k',label='Gas Temp')

    if oplot_t0:
        ax4.plot(atmos.tg1t[:,0]/1000.0, atmos.z1t[:,0]/1e5,color='k',label='Gas Temp (t=0s)',ls='dotted')

  
    srcall = cfline['src']
    srcline = srcall[0].transpose()
            #print(srcline.shape, wl2d.shape, z2d.shape)
            #plvec = np.vectorize(planckfni)
            #b1d = np.zeros((len(atmos.tg1t[:,tind]),1))
            #b1d[:,0] = plvec(np.median(wlprep), atmos.tg1t[:,tind])
    if flip_wl:
        src_flipped = np.flip(srcline[:,:,tind],axis=1)
    else:
        src_flipped = srcline[:,:,tind]
    dnu2dlam = np.zeros_like(src_flipped) + cnst.CCANG / np.median(wl2d)**2
    src_flipped_dlam = src_flipped * dnu2dlam
    src_slice_dlam = src_flipped_dlam[:,15]
    src_trad = np.zeros(len(src_slice_dlam))
    for tt in range(len(src_slice_dlam)):
        src_trad[tt] = trad(lam_rest, src_slice_dlam[tt])
        
    ax4.plot(src_trad/1000, atmos.z1t[:,tind]/1e5,color=bright[4],ls='dashed',label='Source Fctn Temp')

    if source_ha == True:
        srcall_ha = cfline_ha['src']
        srcline_ha = srcall_ha[0].transpose()
        if flip_wl:
            src_flipped_ha = np.flip(srcline_ha[:,:,tind],axis=1)
        else:
            src_flipped_ha = srcline_ha[:,:,tind]
        dnu2dlam_ha = np.zeros_like(src_flipped_ha) + cnst.CCANG / np.median(6562.8)**2
        src_flipped_dlam_ha = src_flipped_ha * dnu2dlam_ha
        src_slice_dlam_ha = src_flipped_dlam_ha[:,25]
        src_trad_ha = np.zeros(len(src_slice_dlam_ha))
        for tt in range(len(src_slice_dlam_ha)):
            src_trad_ha[tt] = trad(6562.8, src_slice_dlam_ha[tt])
        
        ax4.plot(src_trad_ha/1000, atmos.z1t[:,tind]/1e5,color=bright[4],ls='dashed',label=r'H$\alpha$ Source Function')
    
    ax4.set_ylim(ylim)
    ax4.set_xlabel(r'Temp ($10^3$ K)')
    ax4.legend(fontsize=8,frameon=False,loc=(.02,1.03))
    tval = user_time
    plt.title('t = {0:.2f}'.format(user_time)+' s')
    plt.tight_layout()
    if savefig:
        plt.savefig(plot_filename+'.pdf')
        #print('created plot: ',plot_filename+'.pdf')
    #print('The map value range is [',vmin,vmax,'] erg/s/cm3/sr/Ang.')
    #print('Returned height (km), wl (Ang), and contribution function with dimensions:',zprep2.shape, wlprep2.shape, c_I_line.shape)
    tauall = cfline['tau']
    tauline = tauall[0].transpose()
         #   print(tauline.shape, wl2d.shape, z2d.shape)
            #ax1.pcolormesh(wl2d, z2d,np.log10(np.flip(tauline[:,:,tind]/mu,axis=1)),vmax=vmax,vmin=vmin,cmap = xmap,rasterized=True)
            
    ilam_ret = emerg_intline[:,tind]

    if flip_wl:
        tau_return = np.flip(tauline[:,:,tind]/mu,axis=1)
    else:
        tau_return = tauline[:,:,tind]/mu
    ci_dict = {'zpmesh':z2d, 'wlpmesh':wl2d, 'z1d':zprep2, 'wl1d':wlprep2, 'contrib2D':c_I_line[:,:,tind], 'lam0':lam_rest, 'Ilam':ilam_ret, 'tau2D':tau_return} 
    return ci_dict



def color_bright(printc='no'):
    ''' color_bright(printc = no or yes) '''
    bright = [(68,119,170), (102,204,238), (34, 136, 51), (204,187,68), (238,102,119), (170,51,119), (187,187,187)]   
    labels=['blue' ,'cyan', 'green', 'yellow','red','purple', 'light grey']
    for i in range(len(bright)):    
        r, g, b = bright[i]    
        bright[i] = (r / 255., g / 255., b / 255.)
        if printc == 'yes' or printc =='y':
            print(i, labels[i])
    return bright


def make_dataframe(timet, wl, ilam2d, ilamave=np.zeros(1)):
    nt = len(timet)
    nw = len(wl)
    ny = nt*nw
    wl1d = []
    ilam1d = []
    t1d = [] 
    typ = []
    for i in range(nt):
        for j in range(nw):
            t1d.append(timet[i])
            wl1d.append(wl[j])
            ilam1d.append(ilam2d[i, j])
            typ.append('fl')
        for w in range(nw):
            t1d.append(timet[i])
            wl1d.append(wl[w])
            ilam1d.append(ilam2d[-1, w])
            typ.append('post')
        for h in range(nw):
            t1d.append(timet[i])
            wl1d.append(wl[h])
            ilam1d.append(ilam2d[i, h]-ilam2d[-1, h])
            typ.append('bksub')
        if len(ilamave) > 1:
            for w in range(nw):
                t1d.append(timet[i])
                wl1d.append(wl[w])
                ilam1d.append(ilamave[w])
                typ.append('ave')

    d = {'time': t1d, 'wave': wl1d, 'ilam':ilam1d, 'typ':typ}
    datframe = pd.DataFrame(data=d)
    return datframe

def line_movie(df,w1 = 4300, w2=4400, y1=-2500, y2=26000):
    if y2 <= y1:
        y2 = np.max(df["ilam"]) * 1.10
    fig = px.line(df, x="wave", y="ilam",animation_frame="time",\
        color_discrete_sequence=px.colors.qualitative.Set1, \
                  range_x=[w1, w2], color="typ",range_y=[y1,y2],markers=True)
    fig["layout"].pop("updatemenus") # optional, drop animation buttons
    fig.show()
#https://plotly.com/python/sliders/
#   https://plotly.com/python-api-reference/generated/plotly.express.scatter
#   https://plotly.com/python/discrete-color/
#  Set1, Dark2
#fig.add_traces(list(px.scatter(dataframe_ave, x="wave", y="ilam_ave").select_traces()))



def get_lspec(atmos, time=0.0, line = 'Hg', mu = 0.95):
  
    model_dir = np.load('model_dir.tmp.npy')
    model_contrib = np.load('model_contrib.tmp.npy')
    contribf = scipy.io.readsav(str(model_dir)+str(model_contrib),verbose=False,python_dict=True)
    model_file = np.load('model_file.tmp.npy')
    cf = contribf['lcontribf'][0]
    tind = findind(atmos.timet, time)

    if time < 0:
        tind = -99
    
    if line == 'Ca II K' and mu ==0.95:
        cfline = cf.LCAIIKINTT95
        plt_label = 'Ca II K'
        lam_rest = atmos.alamb[17]
        lam0 = vactoair(lam_rest)

    if line == 'Hg' and mu ==0.95:
        cfline = cf.LHGINTT95
        plt_label = 'Hg'
        lam_rest = atmos.alamb[7]
        lam0 = vactoair(lam_rest)
    if line == 'Hg' and mu ==0.77:
        cfline = cf.LHGINTT77
        plt_label = 'Hg'
        lam_rest = atmos.alamb[7]
        lam0 = vactoair(lam_rest)
    if line == 'HeII304' and mu ==0.95:
        cfline = cf.LHEII304INTT95
        plt_label = 'He II 304'
        lam_rest = atmos.alamb[31]
        lam0 = -99 # vactoair(lam_rest)

    if line == 'Ha' and mu ==0.95:
        cfline = cf.LHAINTT95
        plt_label = 'Ha'
        lam_rest = atmos.alamb[2]
        lam0 = vactoair(lam_rest)
    if line == 'Ha' and mu ==0.77:
        cfline = cf.LHAINTT77
        plt_label = 'Ha'
        lam_rest = atmos.alamb[2]
        lam0 = vactoair(lam_rest)
    if line == 'Hb' and mu ==0.95:
        cfline = cf.LHBINTT95
        plt_label = 'Hb'
        lam_rest = atmos.alamb[4]
        lam0 = vactoair(lam_rest)
    if line == 'Hb' and mu ==0.77:
        cfline = cf.LHBINTT77
        plt_label = 'Hb'
        lam_rest = atmos.alamb[4]
        lam0 = vactoair(lam_rest)
    print('Air and vaccuum rest wavelengths = ',lam0,lam_rest)

    contriblineint = cfline['int']
    emerg_intline = contriblineint[0].transpose()
    wl = cfline['lam']
    emerg_intline_wl = wl[0].transpose()

    if tind >= 0:
        ret_wl = emerg_intline_wl[:,tind]
        ret_ilam =  emerg_intline[:,tind]
        print('Returned wl in Ang (air if > 3000Ang) and i_lam[time, wave] in erg/s/cm2/sr/Ang at time [s] = ',atmos.timet[tind])

    else:
        ret_wl = emerg_intline_wl[:,0]
        ret_ilam = np.transpose(emerg_intline)
        print('Returned wl in Ang (air if > 3000Ang) and all i_lam[time, wave] in erg/s/cm2/sr/Ang')
        
    return ret_wl, ret_ilam



def get_lam0_air(atmos, line = 'Ha'):
    if line == 'Ha': 
        lam_rest = atmos.alamb[2]
    if line == 'Hg':
        lam_rest = atmos.alamb[7]
    if line == 'Ca II K':
        lam_rest = atmos.alamb[17]
    if line == 'Ca II 8542':
        lam_rest = atmos.alamb[20]
    if lam_rest > 3000.:
        lam_rest = vactoair(lam_rest)

    return lam_rest



def akcdf_dz(dz,y,norm=False):
    ciprime = np.zeros_like(y)
    zprime = np.zeros_like(y)
    for j in range(len(dz)):
        ciprime[j] = np.sum(y[0:j]*dz[0:j])
        zprime[j] = np.sum(dz[0:j])
    if norm == True:
        ciprime = ciprime / np.sum(y*dz)
    return zprime, ciprime


def radxsnaps():
    # Creates a plot with up to 5 different times and two y-axes.
    return None


def EBapprox(atmos, time=0, hindex1=119, hindex2=128):
    # for H alpha line.

    tind = findind(atmos.timet, time)
    nl_nlLTE1 = atmos.n1t[hindex1,1,0,tind] /atmos.nstart[hindex1,1,0,tind]
    nl_nuLTE1 = atmos.n1t[hindex1,2,0,tind] /atmos.nstart[hindex1,2,0,tind] 
    nl_nu1 = atmos.n1t[hindex1,1,0,tind] /atmos.n1t[hindex1,2,0,tind] 

    nl_nlLTE2 = atmos.n1t[hindex2,1,0,tind] /atmos.nstart[hindex2,1,0,tind]
    nl_nuLTE2 = atmos.n1t[hindex2,2,0,tind] /atmos.nstart[hindex2,2,0,tind] 
    nl_nu2 = atmos.n1t[hindex2,1,0,tind] /atmos.n1t[hindex2,2,0,tind] 

    nlLTE_nuLTE1 = atmos.nstart[hindex1,1,0,tind] /atmos.nstart[hindex1,2,0,tind]
    nlLTE_nuLTE2 = atmos.nstart[hindex2,1,0,tind] /atmos.nstart[hindex2,2,0,tind]

    gu = 2.0 * 3.0**2
    gl = 2.0 * 2.0**2

    B1 = 1./(gu/gl * nlLTE_nuLTE1 - 1.0)
    B2 = 1./(gu/gl * nlLTE_nuLTE2 - 1.0)

    S1 = 1./(gu/gl * nl_nu1 - 1.0)
    S2 = 1./(gu/gl * nl_nu2 - 1.0)
    print(nl_nlLTE1, nl_nuLTE1, nl_nu1)
    print(nl_nlLTE2, nl_nuLTE2, nl_nu2)

    print('Ratio of height index 1 Source to Planck is ',S1/B1)
    print('Ratio of height index 2 Source to Planck is ',S2/B2)
    print('Ratio of Source at height index 1 to height index 2 is ',S1/S2)



    # this is only for H alpha below #
    model_dir = np.load('model_dir.tmp.npy')
    model_contrib = np.load('model_contrib.tmp.npy')
    contribf = scipy.io.readsav(str(model_dir)+str(model_contrib),verbose=False,python_dict=True)
    model_file = np.load('model_file.tmp.npy')
    cf = contribf['lcontribf'][0]
    cfline = cf.LHAINTT95
    contribline = cfline['contribf']
    contriblinewl = cfline['lam']
    
    c_I_line = contribline[0].transpose()  # at mu = 0.95
    c_I_line_wl = contriblinewl[0].transpose()
    src = cfline['src']
    srcline = src[0].transpose()
    srcshape = srcline.shape

    tauall = cfline['tau']
    tauline = tauall[0].transpose()
    shptau = tauline.shape
    src_eb = np.zeros((shptau[1]))
    for tt in range(shptau[1]):
        src_eb[tt] = np.interp(0.95, tauline[:,tt,tind], srcline[:,tt,tind])

    src_eb = np.flip(src_eb)

    print('Ratio of source functions from RADYN: ',srcline[hindex1,0,tind]/srcline[hindex2,0,tind])
    
    return src_eb, c_I_line_wl



def integrate_ci_full(atmos,line='Ha',wave=4200.0,time=0.0,mu=0.95):
    '''  default values for the call are indicated here: 
    your_python_dictionary = adx.ci_image1(atmos,time=18,vmin=-2,vmax=1,xlim=[6560,6566],\
        user_cmap=rnbw_map,savefig=True,user_figsize=(6,6),ci_log=True,user_fontsize=14,\
                                  vel_color='#BBBBBB',oplot_tau1=1)  '''
    model_dir = np.load('model_dir.tmp.npy')
    model_contrib = np.load('model_contrib.tmp.npy')
    contribf = scipy.io.readsav(str(model_dir)+str(model_contrib),verbose=False,python_dict=True)
    model_file = np.load('model_file.tmp.npy')
    cf = contribf['lcontribf'][0]
  #  cf.dtype
    col_rnbw = color_rainbow14()
    if mu < 0.51:
        print('Are you sure you want to look at a mu-value <= 0.5 for a plane parallel flare atmosphere?')
    muind = findind(atmos.zmu, mu)
    user_mu = atmos.zmu[muind]
    if line == 'Ha' and mu==0.95: 
        cfline = cf.LHAINTT95
        plt_label = r'H$\alpha$'
        lam_rest = atmos.alamb[2]
    if line == 'Hg' and mu ==0.95:
        cfline = cf.LHGINTT95
        plt_label=r'H$\gamma$'
        lam_rest = atmos.alamb[7]
        # need to finish
        
    if line == 'Hg' and mu ==0.77:
        cfline = cf.LHGINTT77
        plt_label=r'H$\gamma$'
        lam_rest = atmos.alamb[7]
        # need to finish

    if line == 'Ca II K' and mu ==0.95:
        cfline = cf.LCAIIKINTT95
        plt_label = 'Ca II K'
        lam_rest = atmos.alamb[17]
    if line == 'Ca II 8542' and mu ==0.95:
        cfline = cf.LCAII8542INTT95
        plt_label = 'Ca II 8542'
        lam_rest = atmos.alamb[20]

    if lam_rest > 3000.:
        lam_rest = vactoair(lam_rest)

    contribline = cfline['contribf']
    contriblinewl = cfline['lam']
    c_I_line = contribline[0].transpose()  # at mu = 0.95
    c_I_line_wl = contriblinewl[0].transpose()
    tind = findind(atmos.timet, time)
    print(tind, atmos.timet[tind],atmos.timet[tind+1], lam_rest)
    dzt = atmos.dzt[:,tind]
    nw = len(c_I_line_wl[:,0])
    spec_integral = np.zeros(nw)
    for ww in range(0,nw):
        spec_integral[ww] = np.abs(np.sum(dzt * c_I_line[:,ww,tind]))
    
    return spec_integral

def trad(wl, ilam):
    # can solve for T_rad on your own using the blackbody formula:
    # ilam in erg/s/cm2/sr/ang,  wl in ang.
    inuv_percm = ilam * 1e8
    hh, cc, kb = 6.626e-27, 2.998e10, 1.38e-16
    Trad = (1./(np.log( (2.0 * hh * cc**2) / (inuv_percm * (wl/1e8)**5) + 1.0))) * \
        hh * cc / (kb * wl/1e8)
    return Trad
    
def planckfni(inwave, temp):
    ''' planckfn(wavelength[ang], temp[k])
returns intensity B_lam in units of erg/s/cm2/sr/Angstrom
  can vectorize it along one axis via:
    #  import myfunc as mf
    #  pl = np.vectorize(mf.planckfn)
    #  out = pl(wavearr, singletemperature) '''
    # using my own planck function b/c astropy's keeps changing.
    wl = inwave / 1e8
    Blam = 2.0 * cnst.HPL * (cnst.CC)**2 / (wl)**5 * 1.0 / (np.exp(cnst.CC * cnst.HPL / (wl * cnst.KB * temp)) - 1.0)  * 1e-8  # returns I_lam (intensity per Ang)
    # need to multiply by cnst.PI for flux.
    return Blam # erg/s/cm2/sr/Ang


def planckfni_hz(inwave, temp):
    ''' planckfn(wavelength[ang], temp[k])
returns intensity B_lam in units of erg/s/cm2/sr/Angstrom
  can vectorize it along one axis via:
    #  import myfunc as mf
    #  pl = np.vectorize(mf.planckfn)
    #  out = pl(wavearr, singletemperature) '''
    # using my own planck function b/c astropy's keeps changing.
    wl = inwave / 1e8
    Blam = 2.0 * cnst.HPL * (cnst.CC)**2 / (wl)**5 * 1.0 / (np.exp(cnst.CC * cnst.HPL / (wl * cnst.KB * temp)) - 1.0)  * 1e-8  # returns I_lam (intensity per Ang)
    # need to multiply by cnst.PI for flux.
    nuwave = cnst.CCANG / wl
    Bnu = Blam * cnst.CCANG / (nuwave**2)
    return Bnu # erg/s/cm2/sr/Hz

