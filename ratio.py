from scipy.integrate import odeint, quad
from math import exp, erf, log, pi
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib import ticker 

    # Constants used in model #

"empty array"
V_A_LIST, V_P_LIST = [], []
RATIO_ARRAY = []

for i in range(0, 37):
    print(f"run: {i}")
    F = 96.485333		#pColoumb / fmol
    R = 8.3144598		#pColoumb * mV / (K * fmol)
    T = 310  			#Kelvins
    D =9e-11#5e-10##0.175*5e-10
    dx = 1e-5 #m
    dx_ref = 1e-5 #m
    D2 =0.001*D
    D3 =0.001*D

    A = 2e-10        #m^2 crosssectional area 
    U_KCl = 1.3e-3			#fmol / (s * mV)	Potassium Chloride cotransporter strength
    U_NaKCl2 = 1.3e-3#/2 #35	
    # Parameters used in model

    C_a = 20					#pColoumb / mV 	Membrane Capacitance
    C_p = 20

    P_Na_T =0# 1200e-3

    P_Na_a_T = P_Na_T			#um^3 / s		maximal transient sodium permeability
    P_Na_p_T = P_Na_T

    P_Na_L = 1.51e-3  #(adjust)

    P_Na_a_L = 1.51e-3#*(1+1/3)				#um^3 / s		leak sodium permeability
    P_Na_p_L = 1.51e-3

    P_K_D =0#600e-3

    P_K_a_D = P_K_D				#um^3 / s		maximal delayed rectifier potassium permeability
    P_K_p_D = P_K_D		

    P_K_L = 20.0e-3
    P_K_a_L = P_K_L         	#um^3 / s		leak potassium permeability
    P_K_p_L = P_K_L

    P_Cl_L = 2.5e-3
    P_Cl_a_L = 2.5e-3			#um^3 / s		leak chlroide permeability
    P_Cl_p_L = 2.5e-3
    print(f"ratio: {(P_Na_a_L/P_K_a_L * i/10)/(P_Na_p_L/P_K_p_L)}")
    RATIO_ARRAY.append((P_Na_a_L/P_K_a_L * i/10)/(P_Na_p_L/P_K_p_L))
    "Change"
    # P_K_a_L *= 0.5 #um^3 / s		leak potassium permeability
    # P_K_a_L *= 2 #um^3 / s		leak potassium permeability
    # P_Cl_a_L *= 0.5 #um^3 / s		leak chlroide permeability
    # P_Cl_a_L *= 2 #um^3 / s		leak chlroide permeability
    P_Na_a_L *= i/10 #um^3 / s		leak sodium permeability
    # P_Na_a_L *= 0.5 #um^3 / s		leak sodium permeability

    # P_K_p_L *= 0.5 #um^3 / s		leak potassium permeability
    # P_K_p_L *= 2 #um^3 / s		leak potassium permeability
    # P_Cl_p_L *= 0.5 #um^3 / s		leak chlroide permeability
    # P_Cl_p_L *= 2 #um^3 / s		leak chlroide permeability
    # P_Na_p_L *= 2 #um^3 / s		leak sodium permeability
    # P_Na_p_L *= 0.5 #um^3 / s		leak sodium permeability

    Q_a_pump =54.5#0		#pColoumb / s	Maximal sodium/potassium pump current
    Q_p_pump =54.5	#pColoumb / s	Maximal sodium/potassium pump current


    Na_conce = 152		#mM		Extracellular bath sodium concentration
    K_conce = 3				#mM		Extracellular bath potassium concentration
    Cl_conce = 135.0 # instead of 135.0			#mM		Extracellular bath Chloride concentration
    #cation_conce = 5 #Extracellular cation concentration 

    L_H2O =2#2					#um^3 (s bar)	Effective membrane water permeability

    S_e = 310				#fmol / um^3		Total extracellular solute concentration(original data?)
    V_desired = -70

    # Initial Variables #

    m_a_init = 1.3e-2						#Dimensionless	Initial transient sodium activation gate
    m_p_init = 1.3e-2

    h_a_init = 0.987						#Dimensionless	Initial transient sodium inactavation gate
    h_p_init = 0.987	

    n_a_init = 3e-3						#Dimensionless	Initial delayed rectifier potassium activation gate
    n_p_init = 3e-3		

    W_a_init = 2000.0 
    W_p_init = 2000.0

    dm_anions_init = 0

    Na_a_conci_init = 10.0 					#mM	Initial intracellular sodium concentration
    Na_p_conci_init = 10.0 

    K_a_conci_init = 145.0					#mM	Initial intracellular potassium concentration
    K_p_conci_init = 145.0


    Cl_a_conci_init = 9.75#- 2.0 # 3.8					#mM	Initial intracellular chloride concentration
    Cl_p_conci_init = 9.75


    #N_cat_conci_init = 5 

    N_Na_a_init = W_a_init * Na_a_conci_init	#amol			Initial Intracellular Moles of Sodium
    N_Na_p_init = W_p_init * Na_p_conci_init

    N_K_a_init = W_a_init * K_a_conci_init	#amol			Initial Intracellular Moles of Potassium
    N_K_p_init = W_p_init * K_p_conci_init

    N_Cl_a_init = W_a_init * Cl_a_conci_init	#amol			Initial Intracellular Moles of Chloride
    N_Cl_p_init = W_p_init * Cl_p_conci_init


    A_a_baseline = pi * ((3 * W_a_init) / (4 * pi))**(2/3)						#um^2	Initial Surface Area of Cell
    A_p_baseline = pi * ((3 * W_p_init) / (4 * pi))**(2/3)


    W_a_homeo = 2000.0 
    W_p_homeo = 2000.0 


    #These fonctions might be changed, depend on the shapes of two ends (cylinder /sphere)
    A_a_homeo = pi * ((3 * W_a_homeo) / (4 * pi))**(2/3)	
    A_p_homeo = pi * ((3 * W_p_homeo) / (4 * pi))**(2/3)


    N_a_anions, N_p_anions = 290500, 290500

    xm1, xm2 =0.5, 0.5

    m_a_anions, m_p_anions = xm1 * N_a_anions, xm2 * N_p_anions
    im_a_anions, im_p_anions = (1-xm1) * N_a_anions, (1-xm2) * N_p_anions


    dN_Na_a_init = 0		#fmol	Initial Change In Intracellular Moles of Sodium
    dN_Na_p_init = 0

    dN_K_a_init = 0	        #fmol	Initial Change In Intracellular Moles of Potassium
    dN_K_p_init = 0	


    dN_Cl_a_init = 0		#fmol	Initial Change In Intracellular Moles of Chloride
    dN_Cl_p_init = 0

    N_Na_a_init = 0		#fmol			Initial Intracellular Moles of Sodium
    N_Na_p_init = 0

    N_K_a_init = 0		#fmol			Initial Intracellular Moles of Potassium
    N_K_p_init = 0

    N_Cl_a_init = 0		#fmol			Initial Intracellular Moles of Chloride
    N_Cl_p_init = 0

    # Lists Used For Creating Graphs #

    V_a_TOTAL = []
    V_p_TOTAL = []

    t_a_TOTAL = []
    t_p_TOTAL = []

    A_a_change_TOTAL = []
    A_p_change_TOTAL = []

    I_a_Na_TOTAL = []
    I_p_Na_TOTAL = []

    I_a_K_TOTAL = []
    I_p_K_TOTAL = []

    I_a_Cl_TOTAL = []
    I_p_Cl_TOTAL = []

    I_a_total_TOTAL = []
    I_p_total_TOTAL = []

    h_a_TOTAL = []
    h_p_TOTAL = []

    n_a_TOTAL = []
    n_p_TOTAL = []

    m_a_TOTAL = []
    m_p_TOTAL = []


    # Inner concentrations #

    def Na_a_conc(N_a_Na,W_a):
        return (Na_a_conci_init*W_a_init + N_a_Na)/W_a

    def Na_p_conc(N_p_Na,W_p):
        return (Na_p_conci_init*W_p_init + N_p_Na)/W_p



    def K_a_conc(N_a_K,W_a):
        return (K_a_conci_init*W_a_init + N_a_K)/W_a

    def K_p_conc(N_p_K,W_p):
        return (K_p_conci_init*W_p_init + N_p_K)/W_p




    def Cl_a_conc(N_a_Cl,W_a):
        return (Cl_a_conci_init*W_a_init + N_a_Cl)/W_a

    def Cl_p_conc(N_p_Cl,W_p):
        return (Cl_p_conci_init*W_p_init + N_p_Cl)/W_p



    # Nernst potentials #

    def E_a_Na( N_a_Na,W_a):
        return -R*T/F*log(Na_a_conc(N_a_Na,W_a)/Na_conce) 
    def E_p_Na( N_p_Na,W_p):
        return -R*T/F*log(Na_p_conc(N_p_Na,W_p)/Na_conce) 


    def E_a_K( N_a_K,W_a):
        return -R*T/F*log(K_a_conc(N_a_K,W_a)/K_conce) 
    def E_p_K( N_p_K,W_p):
        return -R*T/F*log(K_p_conc(N_p_K,W_p)/K_conce) 


    def E_a_Cl( N_a_Cl,W_a):
        return R*T/F*log(Cl_a_conc(N_a_Cl,W_a)/Cl_conce) 
    def E_p_Cl( N_p_Cl,W_p):
        return R*T/F*log(Cl_p_conc(N_p_Cl,W_p)/Cl_conce) 
    #######################################

    # Sodium GHK Permeability #

    Max_P_Na_GHK =0# 0.05 # 0.001 # 0.05
    # duration of plateau
    risetime = 0.00043 # 0.00025
    plateau = 100 #0.0001 #0.00025 # 0.0025
    decaytime = 0.001

    def P_Na_GHK_plateau(t):
        if t < risetime:
            P_Na_GHK = Max_P_Na_GHK * t / risetime
        elif t < risetime + plateau: #0.0025:
            P_Na_GHK = Max_P_Na_GHK * 1
        else:
            P_Na_GHK = Max_P_Na_GHK * exp(-(t- risetime - plateau)/decaytime)
        return P_Na_GHK

    # Sodium AChR Permeability (Wang and Rich, 2004)

    AChRx = np.array([0,0.10278,0.12582,0.14886,0.1492,0.1495,0.15086,0.16589,0.168,0.17,0.17189,
            0.18213,0.19365,0.20517,0.21,0.2282,0.25124,0.265,0.27428,0.30755,0.31618,0.33,
            0.35363,0.38818,0.41122,0.42139,0.49654,0.58672,0.66186,0.70695,0.7821,0.90234,
            0.99225,1.06136,1.12919,1.15223,1.21796,1.35322,1.40307,1.53358,1.6,1.67695,1.80411,
            1.83417,1.90932,1.96944,1.98447,2.04459,2.11973,2.19488,2.27003,2.36021,2.42033,2.48044,
            2.58565,2.64577,2.70589,2.81109,2.9163,3.00648,3.0666,3.23192,3.54754,3.6,3.69784,3.75,3.80305,
            3.81808,3.9,3.96837,4.10364,4.28399,4.35914,4.43429,4.50944,4.58459,4.65973,4.73488,4.91524,
            5.06553,5.48636,5.50139,5.63666,5.65169,5.78695,5.80198,6.0, 8.0,  1000.0, 10000000000.0]) #, 1000.0])#/1000.0

    AChRy = np.array([0 ,2.7369,6.59017,10.44344,14,16,20,24.77023,26,26.5,27.29699,31.96053,36.94817,41.43165,
            46.58135,51.06483,55.08015,61.02211,63,66.5,70,73,77,79.98236,81.09873,82.31428,83.16051,
            82.10272,80.41025,78.50621,75.96751,69.8323,65.19949,61.67032,58.14116,55.40426,51.00355,43.81054,
            39.18091,33.86727,30,26.66679,21.17373,19.2697,17.57723,16.51943,15.88475,13.98072,12.28825,
            10.59578,9.53798,8.26863,7.42239,6.99928,5.94148,5.3068,4.67213,4.24901,3.61433,3.40277,3.19121,
            2.7681,2.133,2.05,2,1.92186,1.6,1.49874,1.4,1.28718,1.1,0.86,0.8,0.7,0.65251,0.645,0.64,0.6,0.5,0.4,
            0.3,0.17,0.12,0.08,0.03,0.005,0, 0, 0,0]) #, 0])#/82.31428

    for i in range(len(AChRx)):
        AChRx[i] = AChRx[i]*0.001
        AChRy[i] = AChRy[i]/82.31428
        
    x = AChRx 
    y = AChRy
    f = interp1d(x, y)
    f2 = interp1d(x, y, kind='cubic')


    delay1 = 0.010
    delay = 0.0
    period = 1.0/100.0
    npermax = 600



    def P_Na_GHK(t):
        return Max_P_Na_GHK * f(t)



    def gcat_Hz (t):
        nper = int((t-delay)/period)
        if t < delay:
            return 0.0
        if t > npermax*period:
            return 0.0
        else:      
            return P_Na_GHK(t -delay -nper * period)



    ############################## Currents #######################################

    # conc of mobile anions
    def mAnions_a_conc (W_a, dm_anions):
        if dm_anions > m_a_anions:
            conc = 2*m_a_anions / W_a
        elif dm_anions < - m_a_anions:
            conc = 0.0
        else: 
            conc = (m_a_anions + dm_anions) / W_a       
        return conc
    
    def mAnions_p_conc (W_p, dm_anions):
        if dm_anions > m_p_anions:
            conc = 0.0
        elif dm_anions < - m_p_anions:
            conc = 2*m_p_anions / W_p
        else: 
            conc = (m_p_anions - dm_anions) / W_p       
        return conc

    # Transmembrane Currents #

    def I_a_Na_trans(V_a, m_a, h_a, N_a_Na, N_a_K, N_a_Cl, W_a):

        return I_a_Na_T(V_a, m_a, h_a, N_a_Na, W_a) + I_a_Na_L(V_a, N_a_Na, W_a) \
            + 3 * I_a_pump(N_a_Na, W_a)-F * J_a_NaKCl2(N_a_Na, N_a_K, N_a_Cl, W_a)


    def I_p_Na_trans(V_p, m_p, h_p, N_p_Na, N_p_K, N_p_Cl, W_p, t):
        return I_p_Na_T(V_p, m_p, h_p, N_p_Na, W_p) + I_p_Na_L(V_p, N_p_Na, W_p)\
        + 3 * I_p_pump(N_p_Na, W_p)  + I_Na_p_GHK(V_p, N_p_Na, W_p, t)-F * J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p)




    def I_a_K_trans(V_a, n_a, N_a_Na, N_a_K,  N_a_Cl, W_a):
        return I_a_K_D(V_a, n_a, N_a_K, W_a) + I_a_K_L(V_a, N_a_K, W_a) - 2 * I_a_pump(N_a_Na, W_a) + F * J_a_KCl(N_a_K, N_a_Cl, W_a)-F * J_a_NaKCl2(N_a_Na,N_a_K, N_a_Cl, W_a)


    def I_p_K_trans(V_p, n_p, N_p_Na , N_p_K,  N_p_Cl, W_p, t):
        return I_p_K_D(V_p, n_p, N_p_K, W_p) + I_p_K_L(V_p, N_p_K, W_p) - 2 * I_p_pump(N_p_Na, W_p) + I_K_p_GHK(V_p, N_p_K, W_p, t)- F * J_p_KCl(N_p_K, N_p_Cl, W_p)- F * J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p)


    def I_a_Cl_trans(V_a,N_a_Na,N_a_K, N_a_Cl ,W_a):
        return I_a_Cl_L(V_a, N_a_Cl, W_a) - F * J_a_KCl(N_a_K, N_a_Cl, W_a)- 2 * F * J_a_NaKCl2(N_a_Na,N_a_K, N_a_Cl, W_a)

    def I_p_Cl_trans(V_p,N_p_Na, N_p_K, N_p_Cl, W_p):
        return I_p_Cl_L(V_p, N_p_Cl, W_p)  - F * J_p_KCl(N_p_K, N_p_Cl, W_p)-2* F * J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p)
        #return I_Cl_L(V, N_Cl, W 
        
    R_ext=1e-3#Gigaohm

    def I_a_trans(V_a, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a):
        return I_a_Na_trans(V_a, m_a, h_a, N_a_Na,N_a_K,N_a_Cl,W_a)+I_a_K_trans(V_a, n_a, N_a_Na, N_a_K,  N_a_Cl, W_a)+I_a_Cl_trans(V_a,N_a_Na,N_a_K, N_a_Cl ,W_a)

    def Vext(V_a, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a):
        return -R_ext*I_a_trans(V_a, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a)
    ############Anions flux (diffusion + drift)################3
    def I_a_diff_anions( W_a, W_p, dm_anions):
        return D2*A/dx*(-(-mAnions_p_conc( W_p, dm_anions)+mAnions_a_conc( W_a,dm_anions)))*1e18/(1e3 / F)

    def I_p_diff_anions( W_a, W_p, dm_anions):
        return - I_a_diff_anions( W_a, W_p, dm_anions)
        
    def I_a_drift_anions( V_a, V_p, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a, W_p, dm_anions):
        f = F / (R*T)
        
        V_ext = Vext(V_a, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a)
        
        if V_a - V_p - V_ext > 0:
            
            curr = D2*A/dx*f*(V_a-V_p- V_ext)*(mAnions_p_conc (W_p, dm_anions))*1e18/(1e3 / F)
        else:
            curr = D2*A/dx*f*(V_a-V_p- V_ext)*(mAnions_a_conc (W_a, dm_anions))*1e18/(1e3 / F)
        return curr

    def I_p_drift_anions(V_a, V_p, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a, W_p, dm_anions):
        
        return - I_a_drift_anions( V_a, V_p, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a, W_p, dm_anions)


    def I_a_anions(V_a, V_p, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a, W_p, dm_anions):
        f = F / (R*T)
        
        V_ext = Vext(V_a, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a)
        
        if V_a - V_p - V_ext> 0:
            curr = D2*A/dx*(-(-mAnions_p_conc( W_p, dm_anions)+mAnions_a_conc( W_a,dm_anions))+
                        f*(V_a-V_p- V_ext)*(mAnions_p_conc (W_p, dm_anions)))*1e18/(1e3 / F)
        else:
            curr = D2*A/dx*(-(-mAnions_p_conc( W_p, dm_anions)+mAnions_a_conc( W_a,dm_anions))+
                        f*(V_a-V_p- V_ext)*(mAnions_a_conc (W_a, dm_anions)))*1e18/(1e3 / F)
        return curr

    def I_p_anions( V_a, V_p, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a, W_p, dm_anions):
        return - I_a_anions(V_a, V_p, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a, W_p, dm_anions ) 
    # Sodium Currents #

    def I_a_Na_T(V_a, m_a, h_a, N_a_Na, W_a):			#Transient Sodium Current (based on a model by Kager et al. (2000))
        x_a = (F * V_a) / (R * T)
        return P_Na_a_T * (m_a**3) * h_a * F * x_a * ((Na_a_conc( N_a_Na,W_a) - (Na_conce * exp(-x_a))) / (1 - exp(-x_a)))

    def I_p_Na_T(V_p, m_p, h_p, N_p_Na, W_p):			#Transient Sodium Current (based on a model by Kager et al. (2000))
        x_p = (F * V_p) / (R * T)
        return P_Na_p_T * (m_p**3) * h_p * F * x_p * ((Na_p_conc( N_p_Na,W_p) - (Na_conce * exp(-x_p))) / (1 - exp(-x_p)))





    def I_a_Na_L(V_a, N_a_Na, W_a):					#Sodium Leak Current
        x_a = (F * V_a) / (R * T)
        return P_Na_a_L * F * x_a * ((Na_a_conc(N_a_Na,W_a) - (Na_conce * exp(-x_a))) / (1 - exp(-x_a)))

    def I_p_Na_L(V_p, N_p_Na, W_p):					#Sodium Leak Current
        x_p = (F * V_p) / (R * T)
        return P_Na_p_L * F * x_p * ((Na_p_conc(N_p_Na,W_p) - (Na_conce * exp(-x_p))) / (1 - exp(-x_p)))


    def I_Na_p_GHK(V_p, N_p_Na, W_p, t):			#Sodium GHK Current (ACh Channels)
        x_p = (F * V_p) / (R * T)
        if x_p > 0.001:
            current = gcat_Hz(t) * x_p / (1.0 - exp(-x_p)) * F * (Na_p_conc(N_p_Na,W_p) - Na_conce * exp(-x_p))
        elif x_p < -0.001:
            current = gcat_Hz(t) * x_p / (1.0 - exp(-x_p)) * F * (Na_p_conc(N_p_Na,W_p) - Na_conce * exp(-x_p))
        else:
            current = gcat_Hz(t) *  F * (Na_p_conc(W_p, N_p_Na) - Na_conce * exp(-x_p)) * (1.0 + x_p/2.0 + x_p*x_p/12.0 - x_p**4 /720.0)
        return current

    # Potassium Currents #

    def I_a_K_D(V_a, n_a, N_a_K, W_a):				#Delayed Rectifier Potassium Current (Kager et al. (2000))
        x_a = (F * V_a) / (R * T)
        return P_K_a_D * (n_a**2) * F * x_a * ((K_a_conc(N_a_K,W_a) - (K_conce * exp(-x_a))) / (1 - exp(-x_a)))
        

    def I_p_K_D(V_p, n_p, N_p_K, W_p):				#Delayed Rectifier Potassium Current (Kager et al. (2000))
        x_p = (F * V_p) / (R * T)
        return P_K_p_D * (n_p**2) * F * x_p * ((K_p_conc(N_p_K,W_p) - (K_conce * exp(-x_p))) / (1 - exp(-x_p)))
        

    def I_K_p_GHK(V_p, N_p_K, W_p, t):				#Potassium GHK Current (ACh Channels)
        x_p = (F * V_p) / (R * T)
        if x_p > 0.001:
            current = gcat_Hz(t) * 1.11 * x_p / (1.0 - exp(-x_p)) * F * (K_p_conc( N_p_K,W_p) - K_conce* exp(-x_p))
        elif x_p < -0.001:
            current = gcat_Hz(t) * 1.11 * x_p / (1.0 - exp(-x_p)) * F *(K_p_conc(N_p_K,W_p) - K_conce* exp(-x_p))
        else:
            current = gcat_Hz(t) * 1.11 *  F * (K_p_conc( N_p_K,W_p) - K_conce* exp(-x_p)) * (1.0 + x_p/2.0 + x_p*x_p/12.0 - x_p**4 /720.0)
        return current

    # Diffusion currents
        




    def I_a_K_L(V_a, N_a_K, W_a):					#Potassium Leak Current
        x_a = (F * V_a) / (R * T)
        return P_K_a_L * F * x_a * ((K_a_conc(N_a_K,W_a) - (K_conce * exp(-x_a))) / (1 - exp(-x_a)))

    def I_p_K_L(V_p, N_p_K, W_p):					#Potassium Leak Current
        x_p = (F * V_p) / (R * T)
        return P_K_p_L * F * x_p * ((K_p_conc( N_p_K,W_p) - (K_conce * exp(-x_p))) / (1 - exp(-x_p)))


    # Chloride Currents #



    def I_a_Cl_L(V_a, N_a_Cl, W_a):					#Chloride Leak Current
        x_a = (F * V_a) / (R * T)
        return P_Cl_a_L * F * x_a * ((Cl_a_conc( N_a_Cl,W_a) - (Cl_conce * exp(x_a))) / (1 - exp(x_a)))

    def I_p_Cl_L(V_p, N_p_Cl, W_p):					#Chloride Leak Current
        x_p = (F * V_p) / (R * T)
        return P_Cl_p_L * F * x_p * ((Cl_p_conc( N_p_Cl,W_p) - (Cl_conce * exp(x_p))) / (1 - exp(x_p)))

    # Pump #

    def I_a_pump(N_a_Na, W_a):					#Sodium/Potassium ATP-ase pump (Hamada et al. (2003))
        x1_a = 0.62 / (1 + (6.7 / Na_a_conc( N_a_Na,W_a))**3)
        y1_a = 0.38 / (1 + (67.6 / Na_a_conc( N_a_Na,W_a))**3)
        return Q_a_pump * (x1_a + y1_a)

    def I_p_pump(N_p_Na, W_p):			#Sodium/Potassium ATP-ase pump (Hamada et al. (2003))
        x1_p = 0.62 / (1 + (6.7 / Na_p_conc( N_p_Na,W_p))**3)
        y1_p = 0.38 / (1 + (67.6 / Na_p_conc( N_p_Na,W_p))**3)
        return Q_p_pump * (x1_p + y1_p)



    # ATP Current #
    ATPperpA = 100./1.6/6.022

    def ATP_a_current(N_a_Na, W_a):
        return I_a_pump(N_a_Na, W_a) * ATPperpA   # amol/s

    def ATP_p_current(N_p_Na, W_p):
        return I_p_pump(N_p_Na, W_p) * ATPperpA   # amol/s


    ############################ Diffusion Equation#############


    def I_Na_a_diff(N_a_Na, N_p_Na,W_a,W_p):
        return D*A/dx*(Na_p_conci_init*W_p_init*((1/W_a)-(1/W_p))+((N_a_Na/W_a)-(N_p_Na/W_p)))*1e18*F/1e3 #/10.38

    def I_Na_p_diff(N_a_Na, N_p_Na, W_a, W_p):
        return -I_Na_a_diff(N_a_Na, N_p_Na,W_a,W_p)

    def I_K_a_diff(N_a_K,N_p_K,W_a,W_p):
        return D*A/dx*(K_p_conci_init*W_p_init*((1/W_a)-(1/W_p))+((N_a_K/W_a)-(N_p_K/W_p)))*1e18*F/1e3

    def I_K_p_diff(N_a_K, N_p_K, W_a, W_p):
        return -I_K_a_diff(N_a_K,N_p_K,W_a,W_p)


    def I_Cl_a_diff(N_a_Cl, N_p_Cl,W_a,W_p):
        return - D*A/dx*(Cl_p_conci_init*W_p_init*((1/W_a)-(1/W_p))+((N_a_Cl/W_a)-(N_p_Cl/W_p)))*1e18*F/1e3

    def I_Cl_p_diff(N_a_Cl, N_p_Cl,W_a,W_p):
        return -I_Cl_a_diff(N_a_Cl, N_p_Cl,W_a,W_p)
    ############Drift####################


    def I_Na_a_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_p_Na, N_a_K, N_a_Cl, W_a,W_p):
        f = F / (R*T)
        V_ext = Vext(V_a, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a)
        if V_a - V_p- V_ext> 0:
            drift = D*A*f*(1/dx)*(V_a-V_p-V_ext)*Na_a_conc( N_a_Na,W_a)*1e18/(1e3 / F)
        else:
            drift = D*A*f*(1/dx)*(V_a-V_p-V_ext)*Na_p_conc( N_p_Na,W_p)*1e18/(1e3 / F)
        return drift

    #    return D*A*(F/R*T)*(1/dx)*(V_a - V_p)*(Na_a_conc( N_a_Na,W_a)+Na_p_conc(N_p_Na,W_p))/2*1e18*F/1e3


    def I_Na_p_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_p_Na, N_a_K, N_a_Cl, W_a,W_p):
        return -I_Na_a_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_p_Na, N_a_K, N_a_Cl, W_a,W_p)


    def I_K_a_drift(V_a, V_p,m_a, h_a, n_a,N_a_Na, N_a_K, N_p_K,N_a_Cl,W_a,W_p):
        f = F / (R*T)
        V_ext = Vext(V_a, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a)
        if V_a - V_p -V_ext> 0:
            drift = D*A*f*(1/dx)*(V_a-V_p-V_ext)*K_a_conc( N_a_K,W_a)*1e18/(1e3 / F)
        else:
            drift = D*A*f*(1/dx)*(V_a-V_p-V_ext)*K_p_conc( N_p_K,W_p)*1e18/(1e3 / F)
        return drift
    #    return D*A*(F/R*T)*(1/dx)*(V_a - V_p)* (K_a_conc( N_a_K,W_a)+ K_p_conc( N_p_K,W_p))/2*1e18*F/1e3

    def I_K_p_drift(V_a, V_p,m_a, h_a, n_a,N_a_Na, N_a_K, N_p_K,N_a_Cl,W_a,W_p):
        return -I_K_a_drift(V_a, V_p,m_a, h_a, n_a,N_a_Na, N_a_K, N_p_K,N_a_Cl,W_a,W_p)



    def I_Cl_a_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl,N_p_Cl,W_a,W_p):
        f = F / (R*T)
        V_ext = Vext(V_a, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a)
        if V_a - V_p-V_ext > 0:
            drift = D*A*f*(1/dx)*(V_a-V_p-V_ext)*Cl_p_conc( N_p_Cl,W_p)*1e18/(1e3 / F)
        else:
            drift = D*A*f*(1/dx)*(V_a-V_p-V_ext)*Cl_a_conc( N_a_Cl,W_a)*1e18/(1e3 / F)
        return drift
    #    return -D*A*(F/R*T)*(1/dx)*(V_a - V_p)* (Cl_a_conc( N_a_Cl,W_a)+Cl_p_conc( N_p_Cl,W_p))/2*1e18*F/1e3

    def I_Cl_p_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl,N_p_Cl,W_a,W_p):
        return -I_Cl_a_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl,N_p_Cl,W_a,W_p)
    #################################

    def J_a_KCl(N_a_K, N_a_Cl, W_a):				#Potassium Chloride Flux (Blaesse et al., (2009))
        return 0# (U_KCl * ((R * T) / F) * log((K_a_conc(N_a_K,W_a) * Cl_a_conc(N_a_Cl,W_a)) / (K_conce * Cl_conce)))

    def J_p_KCl(N_p_K, N_p_Cl, W_p):				#Potassium Chloride Flux (Blaesse et al., (2009))
        return 0#(U_KCl * ((R * T) / F) * log((K_p_conc( N_p_K,W_p) * Cl_p_conc( N_p_Cl,W_p)) / (K_conce * Cl_conce)))

    def J_a_NaKCl2(N_a_Na,N_a_K, N_a_Cl, W_a):			
        return 0#(U_NaKCl2 * ((R * T) / F) * log((Na_a_conc(N_a_Na,W_a) * K_a_conc(N_a_K,W_a) * (Cl_a_conc(N_a_Cl,W_a))**2) / (Na_conce * K_conce * ( Cl_conce)**2)))

    def J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p):				
        return 0#(U_NaKCl2 * ((R * T) / F) * log((Na_p_conc( N_p_Na,W_p) * K_p_conc( N_p_K,W_p) * (Cl_p_conc( N_p_Cl,W_p))**2) / (K_conce * (Cl_conce)**2))

    #######Anions concentration#############

    # Total Currents #

    def I_a_Na_tot(V_a,V_p, m_a, h_a,n_a, N_a_Na, N_p_Na,N_a_K, N_a_Cl,W_a, W_p):

        return I_a_Na_T(V_a, m_a, h_a, N_a_Na, W_a) + I_a_Na_L(V_a, N_a_Na, W_a) \
            + 3 * I_a_pump(N_a_Na, W_a) + I_Na_a_diff(N_a_Na, N_p_Na,W_a,W_p) + I_Na_a_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_p_Na, N_a_K, N_a_Cl, W_a,W_p)-F * J_a_NaKCl2(N_a_Na, N_a_K, N_a_Cl, W_a)


    def I_p_Na_tot(V_a,V_p,m_a, m_p, h_a,h_p,n_a, N_a_Na,N_p_Na, N_a_K,N_p_K, N_a_Cl,N_p_Cl,W_a, W_p, t):
        return I_p_Na_T(V_p, m_p, h_p, N_p_Na, W_p) + I_p_Na_L(V_p, N_p_Na, W_p)\
        + 3 * I_p_pump(N_p_Na, W_p) + I_Na_p_diff(N_a_Na, N_p_Na, W_a, W_p) + I_Na_p_GHK(V_p, N_p_Na, W_p, t)+ I_Na_p_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_p_Na, N_a_K, N_a_Cl, W_a,W_p)-F * J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p)




    def I_a_K_tot(V_a,V_p,m_a,h_a ,n_a, N_a_Na, N_a_K, N_p_K, N_a_Cl, W_a, W_p):
        return I_a_K_D(V_a, n_a, N_a_K, W_a) + I_a_K_L(V_a, N_a_K, W_a) - 2 * I_a_pump(N_a_Na, W_a) + I_K_a_diff(N_a_K,N_p_K,W_a,W_p)+I_K_a_drift(V_a, V_p,m_a, h_a, n_a,N_a_Na, N_a_K, N_p_K,N_a_Cl,W_a,W_p)+ F * J_a_KCl(N_a_K, N_a_Cl, W_a)-F * J_a_NaKCl2(N_a_Na, N_a_K, N_a_Cl, W_a)


    def I_p_K_tot(V_a,V_p,m_a,h_a ,n_a, n_p,N_a_Na, N_p_Na, N_a_K, N_p_K, N_a_Cl, N_p_Cl,W_a, W_p, t):
        return I_p_K_D(V_p, n_p, N_p_K, W_p) + I_p_K_L(V_p, N_p_K, W_p) - 2 * I_p_pump(N_p_Na, W_p) + I_K_p_diff(N_a_K, N_p_K, W_a, W_p) + I_K_p_GHK(V_p, N_p_K, W_p, t)+I_K_p_drift(V_a, V_p,m_a, h_a, n_a,N_a_Na, N_a_K, N_p_K,N_a_Cl,W_a,W_p)+ F * J_p_KCl(N_p_K, N_p_Cl, W_p)-F * J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p)


    def I_a_Cl_tot(V_a,V_p,m_a, h_a, n_a, N_a_Na,N_a_K, N_a_Cl, N_p_Cl,W_a,W_p):
        return I_a_Cl_L(V_a, N_a_Cl, W_a) + I_Cl_a_diff(N_a_Cl, N_p_Cl,W_a,W_p)+I_Cl_a_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl,N_p_Cl,W_a,W_p) - F * J_a_KCl(N_a_K, N_a_Cl, W_a)-2*F * J_a_NaKCl2(N_a_Na, N_a_K, N_a_Cl, W_a)

    def I_p_Cl_tot(V_a,V_p,m_a, h_a, n_a, N_a_Na,N_p_Na, N_a_K, N_p_K,N_a_Cl, N_p_Cl, W_a, W_p):
        return I_p_Cl_L(V_p, N_p_Cl, W_p) + I_Cl_p_diff(N_a_Cl, N_p_Cl,W_a,W_p)+I_Cl_p_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl,N_p_Cl,W_a,W_p) - F * J_p_KCl(N_p_K, N_p_Cl, W_p)-2*F * J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p)
        #return I_Cl_L(V, N_Cl, W)
        
        
        


    ########################### Hodgkin-Huxley Gated Rates #########################

    def alpha_m(z):
            x_a = (z + 51.9) / 4
            
            if x_a < -0.001:
                    alpha_a = 1.28e3 * x_a / (1 - exp(-x_a))
            elif x_a > 0.001:
                    alpha_a = 1.28e3 * x_a / (1 - exp(-x_a))
            else:
                    alpha_a = 1.28e3 * (1.0 + x_a /2.0 + x_a*x_a / 12.0 - x_a**4 / 720.0)
            return alpha_a
        

            

    def beta_m(z):
            x_a = (z + 24.91) / 5
            
            if x_a < -0.001:
                    beta_a = 1.4e3 * x_a / (exp(x_a) - 1)
            elif x_a > 0.001:
                    beta_a = 1.4e3 * x_a / (exp(x_a) - 1)
            else:
                    beta_a = 1.4e3 * (1.0 - x_a / 2.0 + x_a*x_a / 12.0 - x_a**4 / 720.0)
            return beta_a


    def alpha_h(z):
        return 0.128e3 * exp(-(z+ 52.92) / (18))




    def beta_h(z):
        return 4e3 / (1 + (exp(-(z + 30) / (5))))


    def alpha_n(z):
            x_a = (z + 34.9) / 5
            
            if x_a < -0.001:
                    alpha_a = 0.08e3 * x_a / (1 - exp(-x_a))
            elif x_a > 0.001:
                    alpha_a = 0.08e3 * x_a / (1 - exp(-x_a))
            else:
                    alpha_a = 0.08e3 * (1.0 + x_a / 2.0 + x_a*x_a / 12.0 - x_a**4 / 720.0)
            return alpha_a


        
    

    def beta_n(z):
        return 0.25e3 * exp(-(z + 50) / (40))



    inh_a_time = 450.
    inh_c_time = 450.
    inh_p_time = 450.

    inh_a_on = 0.002 #150.0
    inh_c_on = 0.002 #150.0
    inh_p_on = 0.002 #150.0

    # Stimulant Current #


    def I_a_stim2(inject_a_strength, inject_a_length, t):
        if t < inh_a_on:
            I_a_s = 0
        elif t < inject_a_length + inh_a_on:
            I_a_s = inject_a_strength
        else:
            I_a_s =0
        return I_a_s
    #    return 0


    def I_p_stim2(inject_p_strength, inject_p_length, t):
        if t < inh_p_on:
            I_p_s = 0
        elif t < inject_p_length + inh_p_on:
            I_p_s = inject_p_strength
        else:
            I_p_s =0
        return I_p_s


    ############# RUN CONDITIONS ################

    """Inject sodium into cell for a short time then return to initial conditions"""


    # Duration of run #

    run_stoptime =10000#10000

    run_no =7

    numpoints =500000# 5000000
    numpoints2 = int(numpoints/10) -1 

    start_time = time.time()

    t = [run_stoptime*float(i)/(numpoints-1) for i in range(numpoints)]
    t2 = [run_stoptime*float(i)/(numpoints-1) for i in range(0,numpoints,10)]


    # Strength of Stim #

    stim_a_length =  run_stoptime  # 0.0005# run_stoptime  #4.0#seconds
    stim_p_length =  run_stoptime

    stim_a_strength = 0.0 # 1000.0 # 88.0 #86.750 # 86.7 # 476.0/2.0 # 500.0 # 1000.0 # 54.750 # 32.4 # 17.9 # 14.45 # 3000.0  # 300.0 # 190.0 # 200.0 # 96.0 # 100.0 # 10.5 # 20.438		#pA
    stim_p_strength = 0.0 



    ############################

    # Differential Equations #



    def vf(k,t):

        m_a, m_p, h_a, h_p, n_a,n_p, N_a_Na, N_p_Na, N_a_K, N_p_K, N_a_Cl, N_p_Cl, W_a, W_p, dm_anions = k

        V_a = (F / C_a) * (N_a_Na + N_a_K - N_a_Cl - dm_anions )
        V_p = (F / C_p) * (N_p_Na + N_p_K - N_p_Cl + dm_anions )
        

        f = [(alpha_m(V_a) * (1 - m_a) - beta_m(V_a) * m_a), # / tau,
            
            (alpha_m(V_p) * (1 - m_p) - beta_m(V_p) * m_p),
            
            (alpha_h(V_a) * (1 - h_a) - beta_h(V_a) * h_a), # / tau
            
            (alpha_h(V_p) * (1 - h_p) - beta_h(V_p) * h_p),

            (alpha_n(V_a) * (1 - n_a) - beta_n(V_a) * n_a), # / tau,
            
            (alpha_n(V_p) * (1 - n_p) - beta_n(V_p) * n_p),
            
            -(1e3 / F) * ( I_a_Na_T(V_a, m_a, h_a, N_a_Na, W_a) + I_a_Na_L(V_a, N_a_Na, W_a) +

                        3 * I_a_pump(N_a_Na, W_a) - I_a_stim2(stim_a_strength, stim_a_length, t)
            
            + I_Na_a_diff(N_a_Na, N_p_Na,W_a,W_p)+ I_Na_a_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_p_Na, N_a_K, N_a_Cl, W_a,W_p) )-1e3 * J_a_NaKCl2(N_a_Na, N_a_K, N_a_Cl, W_a),
            
            -(1e3 / F) * ( I_p_Na_T(V_p, m_p, h_p, N_p_Na, W_p) + I_p_Na_L(V_p, N_p_Na, W_p) +

                        3 * I_p_pump(N_p_Na, W_p) - I_p_stim2(stim_p_strength, stim_p_length, t)
            
            +I_Na_p_diff(N_a_Na, N_p_Na,W_a,W_p) + I_Na_p_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_p_Na, N_a_K, N_a_Cl, W_a,W_p) + I_Na_p_GHK(V_p, N_p_Na, W_p, t))-1e3 * J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p),
            
            -(1e3 / F) * (I_a_K_D(V_a, n_a, N_a_K, W_a) + I_a_K_L(V_a, N_a_K, W_a) - 2 * I_a_pump(N_a_Na, W_a) + 
            I_K_a_diff(N_a_K,N_p_K,W_a,W_p)+I_K_a_drift( V_a, V_p,m_a, h_a, n_a,N_a_Na, N_a_K, N_p_K,N_a_Cl,W_a,W_p)) - 1e3 * J_a_KCl(N_a_K, N_a_Cl, W_a)-1e3* J_a_NaKCl2(N_a_Na, N_a_K, N_a_Cl, W_a),
            
            -(1e3 / F) * (I_p_K_D(V_p, n_p, N_p_K, W_p) + I_p_K_L(V_p, N_p_K, W_p) - 2 * I_p_pump(N_p_Na, W_p) + I_K_p_diff(N_a_K,N_p_K,W_a,W_p) + I_K_p_drift(V_a, V_p,m_a, h_a, n_a,N_a_Na, N_a_K, N_p_K,N_a_Cl,W_a,W_p)+
        
            I_K_p_GHK(V_p, N_p_K, W_p, t)) - 1e3 *J_p_KCl(N_p_K, N_p_Cl, W_p)-1e3* J_p_NaKCl2(N_p_Na,N_p_K, N_p_Cl, W_p),
            
            (1e3 / F) * (I_a_Cl_L(V_a, N_a_Cl, W_a) + I_Cl_a_diff(N_a_Cl, N_p_Cl,W_a,W_p)+I_Cl_a_drift(V_a, V_p,m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl,N_p_Cl,W_a,W_p)) - 1e3 * J_a_KCl(N_a_K, N_a_Cl, W_a)- 1e3 *2* J_a_NaKCl2(N_a_Na, N_a_K, N_a_Cl, W_a),
            
            (1e3 / F) * (I_p_Cl_L(V_p, N_p_Cl, W_p) + I_Cl_p_diff(N_a_Cl, N_p_Cl,W_a,W_p)+I_Cl_p_drift( V_a, V_p,m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl,N_p_Cl,W_a,W_p)) -  1e3 * J_p_KCl(N_p_K, N_p_Cl, W_p)- 1e3 *2* J_p_NaKCl2(N_p_Na, N_p_K, N_p_Cl, W_p),
            
            L_H2O * R * T * (Na_a_conc( N_a_Na,W_a) + K_a_conc( N_a_K,W_a) + Cl_a_conc( N_a_Cl,W_a) + mAnions_a_conc (W_a, dm_anions) + im_a_anions / W_a - S_e),

            L_H2O * R * T * (Na_p_conc( N_p_Na,W_p) + K_p_conc( N_p_K,W_p) + Cl_p_conc( N_p_Cl,W_p) + mAnions_p_conc (W_p, dm_anions) + im_p_anions / W_p - S_e),

            
            (1e3 / F) * I_a_anions(V_a, V_p, m_a, h_a, n_a, N_a_Na, N_a_K, N_a_Cl, W_a, W_p, dm_anions )]

        return f

    # ODE solver parameters #

    abserr = 1.0e-10#  e-12
    relerr = 1.0e-10#  e-12

    k = [ 0.004765939180158293 , 0.004765939180158279 , 0.9961775429327958 , 0.9961775429327958 ,
        0.0011662745253280971 , 0.0011662745253280945 , -4322.3531991890595 , -4322.353199189058 ,
        4295.9062781733455 , 4295.906278173387 , -11.891649659856494 , -11.891649659849646 ,
        1999.8763271913688 , 1999.8763271913692 , -7.833842295287095e-12]

    V_a_init = (F / C_a) * (k[6] + k[8] - k[10] - k[14] )

    V_p_init = (F / C_p) * (k[7] + k[9] - k[11] + k[14] )


        

    Anions_a_conc_init = (1-xm1)*N_a_anions  / k[12]
    Anions_p_conc_init = (1-xm2)*N_p_anions / k[13]
    
    V_ext_init=Vext(V_a_init, k[0], k[2], k[4], k[6], k[8], k[10], k[12])
    ksol = odeint(vf,k,t,atol = abserr, rtol=relerr)
    last = numpoints -1
    ############# PRINT TO SCREEN ################

    V_a = (F / C_a) * (ksol[last,6] + ksol[last,8] - ksol[last,10] - ksol[last,14] )
    V_p = (F / C_p) * (ksol[last,7] + ksol[last,9] - ksol[last,11] + ksol[last,14] )

    ############################ Diffusion Equation#############
    
    def imAnions_a_conc (W_a):
        return(1-xm1)*N_a_anions / W_a
    def imAnions_p_conc (W_p):
        return(1-xm2)*N_p_anions / W_p

    V_ext_end=Vext(V_a, ksol[last,0], ksol[last,2], ksol[last,4], ksol[last,6], ksol[last,8], ksol[last,10], ksol[last,12])

    Na_a_conc_temp = (Na_a_conci_init*W_a_init + ksol[last,6])/ksol[last,12]
    K_a_conc_temp = (K_a_conci_init*W_a_init + ksol[last,8])/ksol[last,12]
    Cl_a_conc_temp = (Cl_a_conci_init*W_a_init + ksol[last,10])/ksol[last,12]
    Na_a_total = Na_a_conci_init*W_a_init + ksol[last,6]
    K_a_total = K_a_conci_init*W_a_init + ksol[last,8]
    Cl_a_total = Cl_a_conci_init*W_a_init + ksol[last,10]
    Anions_a_total = (1-xm1)*(N_a_anions)+ mAnions_a_conc (ksol[last,12], ksol[last,14])*ksol[last,12]
    

    Na_p_conc_temp = (Na_p_conci_init*W_p_init + ksol[last,7])/ksol[last,13]
    K_p_conc_temp = (K_p_conci_init*W_p_init + ksol[last,9])/ksol[last,13]
    Cl_p_conc_temp = (Cl_p_conci_init*W_p_init + ksol[last,11])/ksol[last,13]
    Na_p_total = Na_p_conci_init*W_p_init + ksol[last,7]
    K_p_total = K_p_conci_init*W_p_init + ksol[last,9]
    Cl_p_total = Cl_p_conci_init*W_p_init + ksol[last,11]
    Anions_p_total = (1-xm2)*(N_p_anions )+mAnions_p_conc(ksol[last,13], ksol[last,14])*ksol[last,13]

    
    A_a = pi * ((3 * ksol[last,12]) / (4 * pi))**(2/3)
    A_a_baseline = pi * ((3 * W_a_homeo) / (4 * pi))**(2/3)
    A_a_change = (1 + ((A_a - A_a_baseline) / A_a_baseline)) * 100

    A_p = pi * ((3 * ksol[last,13]) / (4 * pi))**(2/3)
    A_p_baseline = pi * ((3 * W_p_homeo) / (4 * pi))**(2/3)
    A_p_change = (1 + ((A_p - A_p_baseline) / A_p_baseline)) * 100

    #Diffusion
    INa_a_diff = I_Na_a_diff(ksol[last,6],ksol[last,7],ksol[last,12],ksol[last,13]) 

    INa_p_diff = I_Na_p_diff(ksol[last,6],ksol[last,7],ksol[last,12],ksol[last,13]) 


    IK_a_diff = I_K_a_diff(ksol[last,8],ksol[last,9],ksol[last,12],ksol[last,13])

    IK_p_diff = I_K_p_diff(ksol[last,8],ksol[last,9],ksol[last,12],ksol[last,13]) 

    ICl_a_diff = I_Cl_a_diff(ksol[last,10],ksol[last,11],ksol[last,12],ksol[last,13]) 

    ICl_p_diff = I_Cl_p_diff(ksol[last,10],ksol[last,11],ksol[last,12],ksol[last,13]) 



    INa_a_drift = I_Na_a_drift(V_a, V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,7],ksol[last,8],ksol[last,10],ksol[last,12],ksol[last,13])

    INa_p_drift = I_Na_p_drift(V_a, V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,7],ksol[last,8],ksol[last,10],ksol[last,12],ksol[last,13])

    IK_a_drift=I_K_a_drift(V_a, V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,9],ksol[last,10],ksol[last,12],ksol[last,13])

    IK_p_drift=I_K_p_drift(V_a, V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,9],ksol[last,10],ksol[last,12],ksol[last,13])

    ICl_a_drift=I_Cl_a_drift(V_a, V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,11],ksol[last,12],ksol[last,13])

    ICl_p_drift=I_Cl_p_drift(V_a, V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,11],ksol[last,12],ksol[last,13])

    Ia_anions= I_a_anions(V_a,V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12],ksol[last,13] ,ksol[last,14])

    Ia_diff_anions= I_a_diff_anions( ksol[last,12],ksol[last,13] ,ksol[last,14])

    Ia_drift_anions= I_a_drift_anions(V_a,V_p, ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12],ksol[last,13] ,ksol[last,14])

    Ip_anions= I_p_anions( V_a,V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12],ksol[last,13] ,ksol[last,14] )

    Ip_diff_anions= I_p_diff_anions( ksol[last,12],ksol[last,13] ,ksol[last,14])

    Ip_drift_anions= I_p_drift_anions( V_a,V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12],ksol[last,13] ,ksol[last,14])

    Ia_trans=I_a_trans(V_a, ksol[last,0], ksol[last,2], ksol[last,4], ksol[last,6], ksol[last,8], ksol[last,10], ksol[last,12])


    INa_a_T = I_a_Na_T(V_a, ksol[last,0], ksol[last,2], ksol[last,6], ksol[last,12])
    INa_a_L = I_a_Na_L(V_a, ksol[last,6], ksol[last,12])
    INa_a_pump = 3 * I_a_pump(ksol[last,6], ksol[last,12])
    INa_a_NaKCl=F*J_a_NaKCl2(ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12])
    INa_a_trans = I_a_Na_trans(V_a,ksol[last,0], ksol[last,2], ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12])
    INa_a_f = I_a_Na_tot(V_a,V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,7],ksol[last,8],ksol[last,10],ksol[last,12],ksol[last,13])


    IK_a_D = I_a_K_D(V_a, ksol[last,4], ksol[last,8], ksol[last,12])
    IK_a_L = I_a_K_L(V_a, ksol[last,8], ksol[last,12])
    IK_a_pump = - 2.0 * I_a_pump(ksol[last,6], ksol[last,12])
    IK_a_trans = I_a_K_trans(V_a, ksol[last,4], ksol[last,6], ksol[last,8],ksol[last,10],ksol[last,12])
    I_K_a_f = I_a_K_tot(V_a,V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,9],ksol[last,10],ksol[last,12],ksol[last,13])
    IK_a_NaKCl=F*J_a_NaKCl2(ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12])
    IK_a_co = F * J_a_KCl(ksol[last,8],ksol[last,10], ksol[last,12])


    ICl_a_L = I_a_Cl_L(V_a, ksol[last,10], ksol[last,12])
    ICl_a_trans = I_a_Cl_trans( V_a,ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12]) 
    I_Cl_a_f = I_a_Cl_tot( V_a,V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,11],ksol[last,12],ksol[last,13]) 
    ICl_a_co = - F * J_a_KCl(ksol[last,8],ksol[last,10], ksol[last,12])
    ICl_a_NaKCl=2*F*J_a_NaKCl2(ksol[last,6],ksol[last,8],ksol[last,10],ksol[last,12])

    I_a_total = INa_a_f + I_K_a_f +I_Cl_a_f
    
    INa_p_T = I_p_Na_T(V_p, ksol[last,1], ksol[last,3], ksol[last,7], ksol[last,13])
    INa_p_L = I_p_Na_L(V_p, ksol[last,7], ksol[last,13])
    INa_p_pump = 3 * I_p_pump(ksol[last,7], ksol[last,13])
    INa_p_trans = I_p_Na_trans(V_p, ksol[last,1], ksol[last,3], ksol[last,7],ksol[last,9],ksol[last,11],ksol[last,13],t[i])
    INa_p_NaKCl=F*J_p_NaKCl2(ksol[last,7],ksol[last,9],ksol[last,11],ksol[last,13])
    INa_p_f = I_p_Na_tot(V_a,V_p,ksol[last,0],ksol[last,1],ksol[last,2],ksol[last,3],ksol[last,4],ksol[last,6],ksol[last,7],ksol[last,8],ksol[last,9],ksol[last,10],ksol[last,11],ksol[last,12],ksol[last,13],t[i])
    INa_p_GHK = I_Na_p_GHK(V_p, ksol[last,7], ksol[last,13], t[i])


    IK_p_D = I_p_K_D(V_p, ksol[last,5], ksol[last,9], ksol[last,13])
    IK_p_L = I_p_K_L(V_p, ksol[last,9], ksol[last,13])
    IK_p_pump = - 2.0 * I_p_pump(ksol[last,7], ksol[last,13])
    IK_p_trans = I_p_K_trans(V_p,ksol[last,5], ksol[last,7],ksol[last,9],ksol[last,11],ksol[last,13],t[i])
    IK_p_NaKCl=F*J_p_NaKCl2(ksol[last,7],ksol[last,9],ksol[last,11],ksol[last,13])
    I_K_p_f = I_p_K_tot(V_a,V_p,ksol[last,0],ksol[last,2],ksol[last,4],ksol[last,5], ksol[last,6],ksol[last,7], ksol[last,8], ksol[last,9],ksol[last,10],ksol[last,11],ksol[last,12],ksol[last,13],t[i])
    IK_p_GHK = I_K_p_GHK(V_p, ksol[last,9], ksol[last,13], t[i])
    IK_p_co =  F * J_p_KCl(ksol[last,9],ksol[last,11], ksol[last,13])


    ICl_p_L = (V_p, ksol[last,11], ksol[last,13])
    ICl_p_trans = I_p_Cl_trans(V_p,ksol[last,7],ksol[last,9],ksol[last,11],ksol[last,13] )
    ICl_p_NaKCl=2*F*J_p_NaKCl2(ksol[last,7],ksol[last,9],ksol[last,11],ksol[last,13])
    I_Cl_p_f = I_p_Cl_tot(V_a,V_p,ksol[last,0],ksol[last,2],ksol[last,4], ksol[last,6], ksol[last,7],ksol[last,8],ksol[last,9],ksol[last,10], ksol[last,11],ksol[last,12],ksol[last,13] )
    ICl_p_co = - F * J_p_KCl(ksol[last,9],ksol[last,11], ksol[last,13])

    I_p_trans = INa_p_trans + IK_p_trans +ICl_p_trans
    I_p_total = INa_p_f + I_K_p_f +I_Cl_p_f
    I_a_diff_total = INa_a_diff+IK_a_diff+ICl_a_diff+Ia_diff_anions
    I_p_diff_total = INa_p_diff+IK_p_diff+ICl_p_diff+Ip_diff_anions
    I_a_drift_total =INa_a_drift+IK_a_drift+ICl_a_drift+Ia_drift_anions
    I_p_drift_total =INa_p_drift+IK_p_drift+ICl_p_drift+Ip_drift_anions

    Ia_Na_diff_drift_current= INa_a_diff + INa_a_drift 
    Ip_Na_diff_drift_current= INa_p_diff + INa_p_drift


    Ia_K_diff_drift_current= IK_a_diff + IK_a_drift
    Ip_K_diff_drift_current= IK_p_diff + IK_p_drift

    Ia_Cl_diff_drift_current=  ICl_a_diff + ICl_a_drift
    Ip_Cl_diff_drift_current=   ICl_p_diff + ICl_p_drift

    I_Anions_diff_drift_current= Ia_anions + Ip_anions
    
    Icyt = -(Ia_Na_diff_drift_current+Ia_K_diff_drift_current+Ia_Cl_diff_drift_current+Ia_anions)

    print("\n-----------------------------------------------------------\n")
    V_cyt_init, V_cyt_end = V_a_init -V_p_init - V_ext_init, V_a -V_p - V_ext_end
    print(f"V_a: {V_a} mV, V_p: {V_p} mV, V_cyt: {V_cyt_end} mV, V_ext: {V_ext_end} mV")
    V_A_LIST.append(V_a)
    V_P_LIST.append(V_p)

ax = plt.subplot(111)
ax.plot(RATIO_ARRAY, V_A_LIST, color="black")
ax.plot(RATIO_ARRAY, V_P_LIST, color="grey")
ax.spines[['right', 'top']].set_visible(False)
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator()) 
plt.ylabel(r'$V_{cyt}$ (mV)')
plt.xlabel(r'$(P_{Na}:P_K)_a:(P_{Na}:P_K)_p$')
plt.show()
    



    