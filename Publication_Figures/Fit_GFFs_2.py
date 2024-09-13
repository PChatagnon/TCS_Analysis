import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from lmfit import Model
import math 
from matplotlib.colors import LogNorm
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FuncFormatter
import matplotlib.patches as patches
import ROOT


Mv = 3.1
Mn = 0.938
alpha_EM = 1/137.
convert_Gev_to_fm = 0.1973


def Lambda(x, y, z):
    return (x - y - z) * (x - y - z) - 4 * y * z


def T_min(Epho, ma_2=0.0, mb_2 = Mn**2, m1_2=Mv**2, m2_2 = Mn**2):
    s = Mn * Mn + 2 * Mn * Epho
    return ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) - np.sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2)))

def T_max(Epho, ma_2=0.0, mb_2 = Mn**2, m1_2=Mv**2, m2_2 = Mn**2):
    s = Mn * Mn + 2 * Mn * Epho
    return ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) + np.sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2)))

def xi_GPD_approach(Epho, t):
    
    
    W = np.sqrt(Mn**2 + 2*Mn*Epho)
    
    #P_p = np.sqrt( (W**2 - (Mv+Mn)**2) * (W**2 - (Mv-Mn)**2) / (4*W**2))
    ##print(Epho, t, P_p)

    #P = (W**2 - Mn**2)/(2*W)
    #cos_theta = ( t - ( np.sqrt(Mn**2 + P**2) - np.sqrt(Mn**2 + P_p**2) )**2 +P_p**2+P**2) / (2*P_p*P)
    
    #P_plus = (np.sqrt(Mn**2 + P**2) + P)/np.sqrt(2)
    #P_p_plus = (np.sqrt(Mn**2 + P_p**2) + cos_theta*P_p)/np.sqrt(2)
    
    #xi = (P_plus - P_p_plus) / (P_plus + P_p_plus)
    
    xi_1 = (t - Mv**2)/(2*Mn**2+Mv**2 - 2*W**2 -t) # Formula 5.12 in arXiv:hep-ph/0504030v3 with a 1/2 removed to match the formula above
    return xi_1
    


# Define the GPD cross-section
def cross_section_GPD(input_var, C_0, m_A, m_C, debug=False):
    t, Epho = input_var
    
    if(debug):
        print(C_0,m_A,m_C)
    
    A_0 = 0.414
    Conv_factor = 0.3894*(10**6)# GeV-2 to nb
    
    alpha_S = 0.3
    e_q = (2/3)
    Phi_V_2 = 1.0952/(4. * math.pi)
    
    W2 = (Mn**2 + 2*Mn*Epho)
    xi = xi_GPD_approach(Epho, t)
    
    A_t = A_0 / ( 1- (t/(m_A**2)) )**3
    C_t = C_0 / ( 1- (t/(m_C**2)) )**3
    
    if(debug):
        print("A_t, C_t")
        print(A_t, C_t)
    
    H2 = A_t + ((2*xi)**2) * C_t
    E2 = -((2*xi)**2) * C_t
    
    if(debug):
        print("xi, ((2*xi)**2)")
        print(xi, ((2*xi)**2))
        print("H2, E2")
        print(H2, E2)
    
    G_2 = 4.0*( (1-(t/(4*Mn**2)))*E2**2 - 2*E2*(H2+E2) +(1-xi**2)*(H2+E2)**2 ) / (xi**4)
    
    sigma = Conv_factor*( (alpha_EM * e_q**2 ) / ( 4*(W2 - Mn**2)**2 ) ) * ( (16.0*math.pi*alpha_S)**2 / (3*Mv**3) ) * Phi_V_2 * G_2
    
    if(debug):
        print("CS")
        print(sigma)
    
    return sigma


def cross_section_GPD_1(input_var, A_0, C_0, m_A, m_C):
    const1_unitconv= 70.231 #4.6 #
    
    t, Epho = input_var
    s = (Mn**2 + 2*Mn*Epho)
    
    A_t = A_0 / ( 1- (t/(m_A**2)) )**3
    C_t = C_0 / ( 1- (t/(m_C**2)) )**3
    
    xi = xi_GPD_approach(Epho, t)
    
    H2= A_t + 4* xi**2 *C_t
    E2= -xi**2 *C_t*4
    const_2=1/pow(s-Mn*Mn,2);
    myfit =const1_unitconv*const_2*(4 /pow(xi, 4))*((1 - t / (4 * Mn * Mn)) * E2 *E2 - 2 * E2 * (H2 + E2) + (1 - xi * xi) * (H2 + E2) *(H2+E2));
    return myfit

def gpd_based_model(t, s,  A_0, C_0, m_A, m_C):
    const1_unitconv=70.231;
    x = -t
    M=0.938;
    mJpsi=3.097;
    xi = (-x - mJpsi**2)/(2*M**2+mJpsi**2 - 2*s +x);
    A_conf=1.0; #0.8; #5/4
    A = A_0 / ( 1- (t/(m_A**2)) )**3
    C = C_0 / ( 1- (t/(m_C**2)) )**3
    H2= A + 4* xi**2 *C
    E2= -xi**2 *C*4
    const_2=1/pow(s-M*M,2);
    myfit =const1_unitconv*const_2*(4 /pow(xi, 4))*((1 + x / (4 * M * M)) * E2 *E2 - 2 * E2 * (H2 + E2) + (1 - xi * xi) * (H2 + E2) *(H2+E2));
    return myfit


def cross_section_GPD_2(input_var, A_0, C_0, m_A, m_C):
    
    t, Epho = input_var
    
    M=0.938;
    Mj=3.097;
    Q=0.0;
    const1_unitconv=70.231;
    s=M*M-Q*Q+2*M*Epho;
    q_gamma = (1/(2*pow(s,0.5)))*pow(s*s-2*(-Q*Q+M*M)*s+(-Q*Q-M*M)*(-Q*Q-M*M),0.5);
    q_jpsi = (1/(2*pow(s,0.5)))*pow(s*s-2*(Mj*Mj+M*M)*s+(Mj*Mj-M*M)*(Mj*Mj-M*M),0.5);
    E_jpsi = pow(Mj*Mj+q_jpsi*q_jpsi,0.5);
    E_gamma = pow(-Q*Q+q_gamma*q_gamma,0.5);
    P_z = (s-M*M)/(2*pow(s,0.5));
    P_0 = pow(P_z*P_z+M*M,0.5);
    P_plus = (P_0+P_z)/(pow(2,0.5));
    P_prime_z = -(-t+Mj*Mj-2*E_gamma*E_jpsi)/(2*q_gamma);
    K=q_jpsi;
    P_prime_0=pow((-K)*(-K)+M*M,0.5);
    P_prime_plus=(P_prime_0+P_prime_z)/pow(2,0.5);
    xi = (P_plus-P_prime_plus)/(P_plus+P_prime_plus);
    A_conf=1.0; #0.8; #5/4
    # dipole-dipole
    #H2=A0/pow((1+x/(par[0]*par[0])),2)+pow(2*xi,2)*par[1]/pow((1+x/(par[2]*par[2])),2);
    #E2=-pow(2*xi,2)*par[1]/pow((1+x/(par[2]*par[2])),2);
    #tripole-tripole
    H2=A_0/(A_conf*pow((1-t/(m_A*m_A)),3))+pow(2*xi,2)*C_0/(A_conf*pow((1-t/(m_C*m_C)),3));
    E2=-pow(2*xi,2)*C_0/(A_conf*pow((1-t/(m_C*m_C)),3));
    #dipole-tripole
    #H2=A0/pow((1+x/(par[0]*par[0])),3)+pow(2*xi,2)*par[1]/pow((1+x/(par[2]*par[2])),3);
    #E2=-pow(2*xi,2)*par[1]/pow((1+x/(par[2]*par[2])),3);
    const_2=1/pow(s-M*M,2);
    myfit =const1_unitconv*const_2*(4 /pow(xi, 4))*((1 - t / (4 * M * M)) * E2 *E2 - 2 * E2 * (H2 + E2) + (1 - xi * xi) * (H2 + E2) *(H2+E2));
    return myfit

def F(s, t, M, m):
    return (1 / (4096 * M ** 2)) * (
        -9 * M ** 10
        + M ** 8 * (-32 + 68 * m ** 2 + 28 * s + 37 * t)
        + 2
        * M ** 6
        * (256 * m ** 4 + 8 * m ** 2 * (32 * s - 3 * t) + t * (56 - 40 * s - 29 * t))
        + 2
        * M ** 4
        * (
            -136 * m ** 6
            + 64 * s ** 2
            - 56 * s ** 3
            + 8 * m ** 4 * (8 + 27 * s - 64 * t)
            + 3 * t ** 2 * (-24 + 7 * t)
            + 4 * s * t * (-4 + 9 * t)
            - 4 * m ** 2 * (6 * s ** 2 + 32 * s * (1 + 4 * t) + t * (-4 + 25 * t))
        )
        + M ** 2
        * (
            144 * m ** 8
            + 144 * s ** 4
            - 192 * s ** 2 * t
            + 96 * s ** 3 * t
            - 16 * s * (-4 + t) * t ** 2
            + (80 - 13 * t) * t ** 3
            + 96 * m ** 6 * (-6 * s + 7 * t)
            + 32 * m ** 4 * (27 * s ** 2 - 6 * t - 39 * s * t + 8 * t ** 2)
            + 16
            * m ** 2
            * (
                -36 * s ** 3
                + 30 * s ** 2 * t
                + 24 * s * t * (1 + 2 * t)
                + t ** 2 * (-4 + 17 * t)
            )
        )
        - t
        * (2 * m ** 2 - 2 * s - t)
        * (
            64 * m ** 4
            + 8 * m ** 6
            - 8 * s ** 3
            + 76 * m ** 4 * t
            - 16 * t ** 2
            - 90 * m ** 2 * t ** 2
            + t ** 3
            + 4 * s ** 2 * (16 + 6 * m ** 2 + 3 * t)
            - 2 * s * (12 * m ** 4 + 3 * t ** 2 + m ** 2 * (64 + 44 * t))
        )
    )

#Define the holographic QCD cross section
def cross_section_holo(input_var, m_A): # formula 8.4 in PHYS. REV. D 101, 086003 (2020), assuming there is a mistake in the units of N
    t, Epho = input_var
    
    N_norm = 7.768
    Conv_factor = 1.0#0.3894*(10**6) # GeV-2 to nb
    
    
    e = (2/3)
    
    W2 = (Mn**2 + 2*Mn*Epho)
    
    F1 = F(W2, t, Mv, Mn) 
    
    A_0 = 0.414
    A_t = A_0 / ( 1- (t/(m_A**2)) )**3

    
    sigma = Conv_factor * N_norm**2 * ( e**2 / (64*math.pi*(W2-Mn**2)**2) ) * (A_t**2/ (4*Mn**2 * A_0**2)) * F1 * (2*(-t) + 8*Mn**2)
    return sigma

def cross_section_holo_2(input_var, D_0, m_A, m_D, expo=3, A_0=0.414, N=2.032): #Formula in PHYSICAL REVIEW D 106, 086004 (2022)
    t, Epho = input_var
    
    N_e = N #in nb GeV-2
    #N_e = 7.768 * 0.389 * 2/3
       
    s = (Mn**2 + 2*Mn*Epho)
    
    tilde_F = ((s-Mn**2)/2)**4
    
    eta = Mv**2 / ( 2*(s-Mn**2) - Mv**2 + t )

    #A_0 = 0.414
    A_t = A_0 / ( 1- (t/(m_A**2)) )**expo
    D_t = D_0 / ( 1- (t/(m_D**2)) )**expo

    
    sigma = N_e**2 * ( 1. / (64*math.pi*(s-Mn**2)**2) ) * ((A_t+(eta**2)*D_t)**2 / (A_0**2)) * tilde_F * 8
    return sigma

def calculate_max_min_fun(fun, t_values, Epho_value, D_0, m_A, m_D, D_0_err, m_A_err, m_D_err):
    min_values = []
    max_values = []

    for t in t_values:
        results = []
        
        for D_0_val in np.linspace(D_0-D_0_err, D_0+D_0_err, 5):
            for m_A_val in np.linspace(m_A-m_A_err, m_A+m_A_err, 5):
                for m_D_val in np.linspace(m_D-m_D_err, m_D+m_D_err, 5):
                    result = fun((t, Epho_value), D_0=D_0_val, m_A=m_A_val, m_D=m_D_val)
                    results.append(result)
        
        min_values.append(np.min(results))
        max_values.append(np.max(results))

    return np.array(min_values), np.array(max_values)

def cross_section_holo_3(input_var, C_0, m_A, m_C, expo=3, A_0=0.414):
    t, Epho = input_var
    s = (Mn**2 + 2*Mn*Epho)
    const_h = 1 / (64* np.pi * (s - Mn ** 2) ** 2 * 4  * Mn ** 2)
    c2 = 7.768#5.93#7.0 #6.0 #7.768
    # mS=1
    eta = Mv ** 2 / (2 * (s - Mn ** 2) - Mv ** 2 + t)
    # print(eta)
    Ak = A_0 * np.power((1 - t / (m_A ** 2)), -expo)
    Ck= C_0*np.power((1 - t  / (m_C**2)), -expo)
    
    sigma = (
        c2 ** 2 * const_h * 
        (Ak ** 2 + 2 * eta ** 2 * Ak * 4 * Ck + eta ** 4 * (4 * Ck) ** 2) 
        / A_0 ** 2
        * F(s, t, Mv, Mn)
        * (-2 * t + 8 * Mn ** 2)
    )
    
    return sigma

def holographic_model(input_var, m_A, C_0, m_C):
    t, s = input_var
    A0=0.414
    mJpsi = 3.096
    mp=0.938
    M=mp
    x = -t
    const_h = 1 / (64* np.pi * (s - M ** 2) ** 2 * 4  * M ** 2)
    c2 = 7.768#5.93#7.0 #6.0 #7.768
    # mS=1
    eta = mJpsi ** 2 / (2 * (s - M ** 2) - mJpsi ** 2 - x)
    # print(eta)
    Ak = A0 * np.power((1 + x / (m_A ** 2)), -3)
    Ck=C_0*np.power((1 + x / (m_C**2)), -3)
    #Ck = (C_0) * (
    #    (1 + (x / 4) * (1 / m_C ** 2 + 1 / m_A ** 2))
    #    / ((1 + x / m_A ** 2) ** 2 * (1 + x / m_C ** 2) ** 2)
    #)
    # return c2*np.power((1+x/(m*m)),-4)*F(s,-x,mJpsi,M)*(2*x+8*M**2)
    return (
        c2 ** 2
        * const_h
        * (Ak ** 2 + 2 * eta ** 2 * Ak * 4 * Ck + eta ** 4 * (4 * Ck) ** 2)
        / A0 ** 2
        * F(s, -x, mJpsi, mp)
        * (2 * x + 8 * M ** 2)
    )

def Dipole_g(t, A_0, m_A, expo_D):
    return A_0 / ( 1- (t/(m_A**2)) )**expo_D

def Dipole_g_band(t_values, A_0, A_0_error, m_A, m_A_error, expo = 2):
    # Compute function values
    Dipole_g_values = Dipole_g(t_values, A_0, m_A, expo_D = expo)

    # Compute upper and lower bounds for shaded band
    Dipole_g_upper = np.max([Dipole_g(t_values, A_0 + A_0_error, m_A + m_A_error, expo_D = expo),Dipole_g(t_values, A_0 - A_0_error, m_A - m_A_error, expo_D = expo),Dipole_g(t_values, A_0 - A_0_error, m_A + m_A_error, expo_D = expo),Dipole_g(t_values, A_0 + A_0_error, m_A - m_A_error, expo_D = expo)], axis=0)
    Dipole_g_lower = np.min([Dipole_g(t_values, A_0 + A_0_error, m_A + m_A_error, expo_D = expo),Dipole_g(t_values, A_0 - A_0_error, m_A - m_A_error, expo_D = expo),Dipole_g(t_values, A_0 - A_0_error, m_A + m_A_error, expo_D = expo),Dipole_g(t_values, A_0 + A_0_error, m_A - m_A_error, expo_D = expo)], axis=0)#Dipole_g(t_values, A_0 - A_0_error, m_A - m_A_error)
    
    return Dipole_g_values, Dipole_g_upper, Dipole_g_lower

def model_dipole(t_values, sigma_0, m_S):
    return Dipole_g(t_values, sigma_0, m_S, 4)


################################################
################ Integrated CS #################
################################################
def Plot_int_CS():
    with open("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_Latex_Table.txt", "r") as file:
        lines = file.readlines()
        data_HB_int = np.array([[float(value) for value in line.split("&")] for line in lines])

    systematics_CLAS12 = np.array([
        [130.453, 16.9437, 2.9274, 1.55966, 6.02919, 7.81237, 13.7265, 33.8826, 19.5421, 32.7588],
        [ 0, 41.6628, 148.573, 54.8454, 34.6713, 134.213, 16.218, 31.8332, 36.9747, 88.9859],
    ])

    data_HD_int = np.array([
        [ 8.29,  8.47,  8.65,  8.83,  9.01,  9.19,  9.37,  9.55,
             9.73,  9.91, 10.09, 10.27, 10.45, 10.63, 10.81, 10.99,
            11.17, 11.35],
        [ 0.043,  0.136,  0.249,  0.326,  0.206,  0.2  ,  0.489,  0.71 ,
           0.507,  0.683,  0.829,  0.848,  1.321,  0.981,  1.151,  1.114,
           1.594,  1.791],
        [ 0.012,  0.022,  0.029,  0.048,  0.059,  0.06 ,  0.087,  0.134,
           0.08 ,  0.1  ,  0.119,  0.123,  0.193,  0.134,  0.14 ,  0.126,
           0.208,  0.344],
        [ 0.027,  0.026,  0.029,  0.016,  0.056,  0.018,  0.019,  0.064,
           0.019,  0.116,  0.064,  0.059,  0.067,  0.104,  0.051,  0.034,
           0.144,  0.026]
    ])


    fig = plt.figure(figsize=(10, 7))

    points_HB = plt.errorbar(data_HB_int[0],data_HB_int[2], xerr=data_HB_int[1], yerr=data_HB_int[3], fmt='o', capsize=0, color='#00b4d8', label='CLAS12 (this work)')
    points_HD = plt.errorbar(data_HD_int[0],data_HD_int[1], yerr=data_HD_int[2], fmt='o', capsize=0, color='#ce0037', label='GlueX (2023)')

    plt.xlabel(r'$E_{\gamma}~[GeV]$')
    plt.ylabel(r'$\sigma~[nb]$ ')
    plt.gca().set_ylim(top=10)
    plt.gca().set_ylim(bottom=0.01)
    plt.gca().set_xlim(right=11.5)


    for x, y, yerr, yfill in zip(data_HD_int[0], data_HD_int[1], data_HD_int[2], data_HD_int[3]):
        plt.fill_between([x - 0.03, x + 0.03], y - np.sqrt(yfill**2 + (y*0.2)**2), y + np.sqrt(yfill**2 + (y*0.2)**2), color='#ce0050', alpha=0.3)

    for x, y, syst_low, syst_high in zip(data_HB_int[0], data_HB_int[2], systematics_CLAS12[0], systematics_CLAS12[1]):
        #plt.fill_between([x - 0.03, x + 0.03], y*(1-0.2), y *(1+0.2), color='#00b4d0', alpha=0.3)
        plt.fill_between([x - 0.03, x + 0.03], y*(1-(syst_low/100.)), y *(1+(syst_high/100.)), color='#00b4d0', alpha=0.3)

    plt.legend(loc='upper left', ncol=2,  frameon=False, prop={'size': 15})
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().set_yscale('log')
    #plt.show()
    fig.savefig("Integrated_CS.png")

Plot_int_CS()

################################################
################ CLAS12 ########################
################################################
with open("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_8.20_9.28_Latex_Table.txt", "r") as file:
    lines = file.readlines()
    data_1 = np.array([[float(value) for value in line.split("&")] for line in lines])
     
with open("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_9.28_10.00_Latex_Table.txt", "r") as file:
    lines = file.readlines()
    data_2 = np.array([[float(value) for value in line.split("&")] for line in lines]) 
    
with open("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_10.00_10.60_Latex_Table.txt", "r") as file:
    lines = file.readlines()
    data_3 = np.array([[float(value) for value in line.split("&")] for line in lines]) 

# Concatenate data from both files
t_HB = np.concatenate((data_1[0], data_2[0], data_3[0]))
E_gamma_HB = np.concatenate((data_1[2], data_2[2], data_3[2]))
dsdt_nb_HB = np.concatenate((data_1[4], data_2[4], data_3[4]))
err_dsdt_nb_HB = np.concatenate((data_1[5], data_2[5], data_3[5]))

color_CLAS12 = '#00b4d8'#'ffb81c''#ffb81c'
label_HallB = 'CLAS12'
################################################  

################################################
################ HallC  ########################
################################################
data = np.genfromtxt("Results_HallC/41586_2023_5730_MOESM2_ESM(1)/Results_HallC.csv", delimiter=",", skip_header=1, dtype=float)

E_gamma_HC = data[:, 8] 
t_HC = data[:, 10]
dsdt_nb_HC = data[:, 11]
err_dsdt_nb_HC = data[:, 14]
color_hall_c = '#154734'
label_HallC = 'Hall C - 007'
################################################

################################################
################ Hall D  #######################
################################################
data_HD = np.array([
    [0.77, 1.00, 0.92, 9.14, 0.313, 0.092, 0.120],
    [1.00, 1.50, 1.25, 8.96, 0.170, 0.018, 0.008],
    [1.50, 2.00, 1.72, 8.80, 0.097, 0.010, 0.040],
    [2.00, 2.50, 2.24, 8.77, 0.045, 0.007, 0.003],
    [2.50, 3.50, 2.94, 8.78, 0.018, 0.003, 0.009],
    [3.50, 4.50, 3.92, 8.95, 0.030, 0.006, 0.004],
    [4.50, 5.75, 4.95, 9.10, 0.033, 0.013, 0.012]
])

data_HD_1 = np.array([
    [0.49, 0.77, 0.69, 10.00, 0.813, 0.088, 0.092],
    [0.77, 1.00, 0.87, 9.85, 0.499, 0.061, 0.016],
    [1.00, 1.50, 1.21, 9.83, 0.401, 0.037, 0.010],
    [1.50, 2.00, 1.71, 9.83, 0.231, 0.027, 0.006],
    [2.00, 2.50, 2.24, 9.82, 0.120, 0.021, 0.007],
    [2.50, 3.50, 2.97, 9.84, 0.075, 0.011, 0.005],
    [3.50, 4.50, 3.89, 9.86, 0.026, 0.008, 0.006],
    [4.50, 5.75, 5.06, 9.76, 0.019, 0.005, 0.002],
    [5.75, 8.10, 6.37, 9.93, 0.009, 0.004, 0.003]
])

data_HD_2 = np.array([
    [0.35, 0.49, 0.46, 10.96, 1.611, 0.187, 0.139],
    [0.49, 0.77, 0.60, 10.87, 1.150, 0.084, 0.109],
    [0.77, 1.00, 0.88, 10.85, 1.015, 0.089, 0.023],
    [1.00, 1.50, 1.18, 10.86, 0.529, 0.042, 0.023],
    [1.50, 2.00, 1.69, 10.86, 0.242, 0.029, 0.008],
    [2.00, 2.50, 2.24, 10.83, 0.170, 0.025, 0.003],
    [2.50, 3.50, 2.87, 10.82, 0.072, 0.012, 0.008],
    [3.50, 4.50, 3.92, 10.81, 0.051, 0.009, 0.002],
    [4.50, 5.75, 4.93, 10.78, 0.016, 0.005, 0.001],
    [5.75, 8.10, 6.97, 10.70, 0.0058, 0.0026, 0.0008],
    [8.10, 10.30, 8.36, 10.70, 0.0047, 0.0024, 0.0002]
])

t_HD = np.concatenate((data_HD[:, 2], data_HD_1[:, 2], data_HD_2[:, 2]))

E_gamma_HD = np.concatenate((data_HD[:, 3], data_HD_1[:, 3], data_HD_2[:, 3]))
dsdt_nb_HD = np.concatenate((data_HD[:, 4], data_HD_1[:, 4], data_HD_2[:, 4]))
err_dsdt_nb_HD = np.sqrt(np.concatenate((data_HD[:, 5], data_HD_1[:, 5], data_HD_2[:, 5]))**2 + np.concatenate((data_HD[:, 6], data_HD_1[:, 6], data_HD_2[:, 6]))**2)

color_glueX = '#ce0037'
label_HallD = 'GlueX'
################################################

################################################
################ Do the fit  ###################
################################################


t_data_B = -t_HB #np.concatenate((-t_HC, -t_HB, -t_HD))
E_data_B = E_gamma_HB #np.concatenate((E_gamma_HC, E_gamma_HB, E_gamma_HD))
cs_data_B = dsdt_nb_HB #np.concatenate((dsdt_nb_HC, dsdt_nb_HB, dsdt_nb_HD))
sigma_cs_B = err_dsdt_nb_HB #np.concatenate((err_dsdt_nb_HC, err_dsdt_nb_HB, err_dsdt_nb_HD))

t_data_C = -t_HC #np.concatenate((-t_HC, -t_HB, -t_HD))
E_data_C = E_gamma_HC #np.concatenate((E_gamma_HC, E_gamma_HB, E_gamma_HD))
cs_data_C = dsdt_nb_HC #np.concatenate((dsdt_nb_HC, dsdt_nb_HB, dsdt_nb_HD))
sigma_cs_C = err_dsdt_nb_HC #np.concatenate((err_dsdt_nb_HC, err_dsdt_nb_HB, err_dsdt_nb_HD))

t_data_D = -t_HD #np.concatenate((-t_HC, -t_HB, -t_HD))
E_data_D = E_gamma_HD #np.concatenate((E_gamma_HC, E_gamma_HB, E_gamma_HD))
cs_data_D = dsdt_nb_HD #np.concatenate((dsdt_nb_HC, dsdt_nb_HB, dsdt_nb_HD))
sigma_cs_D = err_dsdt_nb_HD #np.concatenate((err_dsdt_nb_HC, err_dsdt_nb_HB, err_dsdt_nb_HD))

t_data_D_C = np.concatenate((-t_HC, -t_HD))
E_data_D_C = np.concatenate((E_gamma_HC, E_gamma_HD))
cs_data_D_C = np.concatenate((dsdt_nb_HC, dsdt_nb_HD))
sigma_cs_D_C = np.concatenate((err_dsdt_nb_HC, err_dsdt_nb_HD))

t_data = np.concatenate((-t_HC, -t_HB, -t_HD))
E_data = np.concatenate((E_gamma_HC, E_gamma_HB, E_gamma_HD))
cs_data = np.concatenate((dsdt_nb_HC, dsdt_nb_HB, dsdt_nb_HD))
sigma_cs = np.concatenate((err_dsdt_nb_HC, err_dsdt_nb_HB, err_dsdt_nb_HD))



experiment_array = np.concatenate((np.full_like(E_gamma_HC, 1), np.full_like(E_gamma_HB, 2),np.full_like(E_gamma_HD, 3)))


model_GPD = Model(cross_section_GPD)
params_GPD = model_GPD.make_params(C_0 = -10000.2, m_A = 2.0, m_C = 1.25, debug=False)
params_GPD['m_C'].min = 0
params_GPD['m_A'].min = 0
params_GPD['debug'].vary = False


model_holo = Model(cross_section_holo_2)
params_holo = model_holo.make_params(D_0=-1.8, m_A=1.612, m_D=0.963)
params_holo['m_A'].min = 0
params_holo['D_0'].min = -3
#params_holo['D_0'].vary = False
params_holo['A_0'].vary = False
params_holo['N'].vary = False
params_holo['expo'].vary = False


# Perform the fit with error bars
#result_GPD = model_GPD.fit(cs_data, input_var = (t_data, E_data), params=params_GPD, weights=1/sigma_cs)
#print(result_GPD.fit_report())
print("///////////////////")
print("Results all data")
result_holo = model_holo.fit(cs_data, input_var = (t_data, E_data), params=params_holo, weights=1/sigma_cs)
print(result_holo.fit_report())
print("///////////////////\n")

print("///////////////////")
print("Results C")
result_holo_C = model_holo.fit(cs_data_C, input_var = (t_data_C, E_data_C), params=params_holo, weights=1/sigma_cs_C)
print(result_holo_C.fit_report())
print("///////////////////\n")

print("///////////////////")
print("Results C and D")
result_holo_D_C = model_holo.fit(cs_data_D_C, input_var = (t_data_D_C, E_data_D_C), params=params_holo, weights=1/sigma_cs_D_C)
print(result_holo_D_C.fit_report())
print("///////////////////\n")

print("Results B")
params_holo['D_0'].value = -1.01331202
params_holo['D_0'].vary = False
result_holo_B = model_holo.fit(cs_data_B, input_var = (t_data_B, E_data_B), params=params_holo, weights=1/sigma_cs_B)
print(result_holo_B.fit_report())
print("///////////////////\n")

print("Results D")
params_holo['D_0'].vary = True
result_holo_D = model_holo.fit(cs_data_D, input_var = (t_data_D, E_data_D), params=params_holo, weights=1/sigma_cs_D)
print(result_holo_D.fit_report())
print("///////////////////\n")



################################################
### Do the 1D fits #############################
################################################
model_dipole_fit = Model(model_dipole)
params_dipole = model_dipole_fit.make_params(sigma_0=2.0, m_S=1.87)

def fit_1D(mask_hallB, E_1, E_2, debug=False):

    print(params_dipole)
    print(params_GPD)
    mask = (E_data >= E_1) & (E_data <= E_2)
    result_dipole = model_dipole_fit.fit(cs_data[mask][mask_hallB], t_values=t_data[mask][mask_hallB], params=params_dipole, weights=1/sigma_cs[mask][mask_hallB])
    
    if debug:
        print(result_dipole.fit_report())
        
    return result_dipole

def dipole_model_ROOT(x, params):
    return params[0] / (1 - (x[0]/(params[1]**2)))**4


def fit_1D_Root(mask_hallB, E_1, E_2, debug=False):

    # Parameters for the dipole model
    params_dipole = [2.0, 1.50]

    # Define the fit function using ROOT's TF1
    fit_function = ROOT.TF1("fit_function", dipole_model_ROOT, E_1, E_2, len(params_dipole))
    fit_function.SetParameter(0, params_dipole[0])
    fit_function.SetParameter(1, params_dipole[1])


    # Create a TGraphErrors object with the data to be fitted
    graph = ROOT.TGraphErrors()
    point_index = 0
    mask = (E_data >= E_1) & (E_data <= E_2)
    print(E_1,E_2)
    for e, t, cs, sigma in zip(E_data[mask][mask_hallB], t_data[mask][mask_hallB], cs_data[mask][mask_hallB], sigma_cs[mask][mask_hallB]):
        graph.SetPoint(point_index, t, cs)
        graph.SetPointError(point_index, 0, sigma)
        point_index += 1
    
    print(point_index)
    # Draw the histogram and the fit function
    canvas = ROOT.TCanvas("canvas", "Fit Example", 800, 600)
    graph.Fit(fit_function)
    graph.Draw("AP")
    canvas.SaveAs("fit_example"+ f"{E_1}"+".png")

    # Perform the fit
    

    if debug:
        fit_function.Print("V")

    return fit_function

################################################
#########   Define the 1D fit table herevv #####
################################################
Eg_values_CLAS12 = np.array([])
eg_errors_up_CLAS12 = np.array([])
eg_errors_down_CLAS12 = np.array([])
mA_values_CLAS12 = np.array([])
mA_errors_CLAS12 = np.array([])
mass_radius_values_CLAS12 = np.array([])
mass_radius_errors_CLAS12 = np.array([])



################################################
################   3D plot   ###################
################################################
# Generate meshgrid for 3D plot
t_min, t_max = -8 , -0.5#np.min(x_data), np.max(x_data)
Epho_min, Epho_max = 8.2, 11#np.min(y_data), np.max(y_data)
t_range = np.linspace(t_min, t_max, 100)
Epho_range = np.linspace(Epho_min, Epho_max, 100)
t_plot, Epho_plot = np.meshgrid(t_range, Epho_range)

cs_from_holo_model = cross_section_holo_2((t_plot, Epho_plot), result_holo.params['D_0'],result_holo.params['m_A'],result_holo.params['m_D'])
#cs_from_holo_model = cross_section_holo_2((t_plot, Epho_plot), D_0=-1.93, m_A=1.641, m_D=1.7, A_0=0.429, expo=3)

#cs = cross_section_GPD((t_plot, Epho_plot), -0.45*1000, 2000, 1.0)

# Create 3D plot
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.errorbar(t_HC, E_gamma_HC, dsdt_nb_HC, zerr=err_dsdt_nb_HC, fmt='o', color=color_hall_c, label=label_HallC)
ax.errorbar(t_HD, E_gamma_HD, dsdt_nb_HD, zerr=err_dsdt_nb_HD, fmt='o', color=color_glueX, label=label_HallD)
ax.errorbar(t_HB, E_gamma_HB, dsdt_nb_HB, zerr=err_dsdt_nb_HB, fmt='o', color=color_CLAS12, label=label_HallB)
surf = ax.plot_surface(-t_plot, Epho_plot, cs_from_holo_model, cmap='viridis', alpha=0.5, label='Fit')
#surf._edgecolors2d = surf._edgecolor3d
#surf._facecolors2d = surf._facecolor3d

ax.set_xlabel(r'$-t~[GeV^2]$ ')
ax.set_ylabel(r'$E_{\gamma}~[GeV]$')
ax.set_zlabel(r'$\frac{d\sigma}{dt}~[nb/GeV^2]$')
#ax.set_zscale('log')
#ax.set_zlim(bottom=7e-4)
#ax.set_zlim(top=2)
#ax.legend()

# Save the plot as PNG
plt.savefig("3d_plot_with_errorbars_and_fit_lmfit.png")

# Show the plot
#plt.show()

################################################
##############   Slice plot   ##################
################################################
def Plot_Slice(name = "Slice_fit.png", nb_columns = 3, bin_limits = [8.2, 9.10, 9.25, 9.40, 9.55, 9.70, 9.85, 10.0, 10.15, 10.30, 10.45, 10.6, 11.5], log = True, only_CLAS12=False, only_Model = True, only_Fit = False, fit_CLAS12_1D = False):
    t_data = np.concatenate((-t_HC, -t_HB, -t_HD))
    E_data = np.concatenate((E_gamma_HC, E_gamma_HB, E_gamma_HD))
    cs_data = np.concatenate((dsdt_nb_HC, dsdt_nb_HB, dsdt_nb_HD))
    sigma_cs = np.concatenate((err_dsdt_nb_HC, err_dsdt_nb_HB, err_dsdt_nb_HD))

    #bin_limits = [8.2, 9.10, 9.25, 9.40, 9.55, 9.70, 9.85, 10.0, 10.15, 10.30, 10.45, 10.6, 11.5]
    #nb_columns = 3
    fig, axs = plt.subplots(math.ceil((len(bin_limits)-1)/nb_columns), nb_columns, figsize=(3.75*nb_columns, (math.ceil((len(bin_limits)-1)/nb_columns))*3)) #math.ceil((len(bin_limits)-1)/2)

    for i in range(len(bin_limits)-1):

        row = i // nb_columns
        col = i % nb_columns

        mask = (E_data >= bin_limits[i]) & (E_data <= bin_limits[i+1])
        t_display = t_data[mask]
        E_display = E_data[mask]
        cs_display = cs_data[mask]
        sigma_cs_display = sigma_cs[mask]
        experiment_array_display = experiment_array[mask]

        mask_hallC = experiment_array_display==1
        mask_hallB = experiment_array_display==2
        mask_hallD = experiment_array_display==3

        tmin=0
        if(not log):
            tmin=0.1
        tmax = math.ceil(max(-t_display)+0.2)
        
        axs[row, col].set_xlim(tmin, tmax)
        axs[row, col].xaxis.set_major_locator(MultipleLocator(1))
        
        if(log):
            axs[row, col].set_yscale('log')
            axs[row, col].set_ylim(bottom=7e-4)
            axs[row, col].set_ylim(top=2)
        else:
            axs[row, col].set_yscale('linear')
            axs[row, col].set_ylim(bottom=0)
            axs[row, col].set_ylim(top=math.ceil(max(cs_display)*1.1))
            
            
        if not only_CLAS12:
            axs[row, col].errorbar(-t_display[mask_hallC], cs_display[mask_hallC], yerr=sigma_cs_display[mask_hallC], fmt='o', capsize=0, color = color_hall_c, label=label_HallC)
            axs[row, col].errorbar(-t_display[mask_hallD], cs_display[mask_hallD], yerr=sigma_cs_display[mask_hallD], fmt='o', capsize=0, color = color_glueX, label=label_HallD)

        axs[row, col].errorbar(-t_display[mask_hallB], cs_display[mask_hallB], yerr=sigma_cs_display[mask_hallB], fmt='o', capsize=0, color = color_CLAS12, label=label_HallB)
        #axs[row, col].set_title(f"E photon in [{bin_limits[i]}, {bin_limits[i+1]}]")
        axs[row, col].annotate(r'$E_{\gamma}$ in ' + f"[{bin_limits[i]:.2f}, {bin_limits[i+1]:.2f}] GeV", xy=(1, 1), xycoords='axes fraction', ha='right', va='top')
        axs[row, col].set_ylabel(' ')

       
            
        t_values = np.linspace(tmin,  8.0+0.2, 100)
        Epho_value = np.mean(E_display) #(bin_limits[i] + (bin_limits[i+1]-bin_limits[i])/2)
        s_values = np.full_like(t_values, Mn**2 + 2*Mn*Epho_value)
        

        #axs[row, col].plot(t_values, cross_section_GPD((-t_values, Epho_value), C_0 = -0.48, m_A=1.64, m_C=1.07),label='GPD N=0.1') #1000*4
        
        if only_Model:
            axs[row, col].plot(t_values, cross_section_holo((-t_values, Epho_value), m_A=1.5),label='Holo') #*1000*4

            axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.275, m_A=1.612, m_D=0.963),label='Holo 2') #*1000*4
            axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.93, m_A=1.641, m_D=1.7, A_0=0.429, expo=3, N=2.549),label='Tripole Lattice') #*1000*4
            axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-10, m_A=1.13, m_D=0.48, A_0=0.58, expo=2, N=3.244),label='Dipole Lattice') #*1000*4

            axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.8, m_A=1.575, m_D=1.12),label='HC Holo fit 2') #*1000*4
            axs[row, col].plot(t_values, holographic_model((-t_values, s_values), m_A=1.575, C_0=-.45, m_C=1.12),label='HC Holo fit') #*1000*4 
            axs[row, col].plot(t_values, cross_section_GPD_1((-t_values, Epho_value), C_0=-0.45, m_A=1.575, m_C=1.12, A_0=0.414),label='HC GPD - Holo fit') #*1000*4
            axs[row, col].plot(t_values, cross_section_GPD_1((-t_values, Epho_value), C_0=-0.2, m_A=1.7, m_C=1.28, A_0=0.414),label='HC GPD - GPD fit ') #*1000*4
        
        if only_Fit:
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=result_holo.params['D_0'], m_A=result_holo.params['m_A'], m_D=result_holo.params['m_D']),label='Fit Hall B') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.19, m_A=1.87, m_D=1.56),label=' My HC Holo fit') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.8, m_A=1.55, m_D=0.63),label=' My HB Holo fit with fixed D') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-3.0, m_A=1.55, m_D=0.57),label=' My HB Holo fit with fixed D') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-0.7, m_A=1.78, m_D=1.54),label=' My HB Holo fit without bin 3') #*1000*4

            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.05, m_A=2.04, m_D=1.96),label=' My HD Holo fit') #*1000*4

            axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.0, m_A=1.86, m_D=1.66), color='green', label='Combined fit') #*1000*4
            min_fit, max_fit = calculate_max_min_fun(cross_section_holo_2, -t_values, Epho_value, D_0=-1.0, m_A=1.86, m_D=1.66, D_0_err=0.101, m_A_err=0.057, m_D_err=0.176)
            axs[row, col].fill_between(t_values, min_fit, max_fit, color='green', alpha=0.5)
            
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.0, m_A=1.55, m_D=0.71),label='CLAS12 only with fixed D(0)') #*1000*4
            #min_fit_1, max_fit_2 = calculate_max_min_fun(cross_section_holo_2, -t_values, Epho_value, D_0=-1.0, m_A=1.55, m_D=0.71, D_0_err=0.101, m_A_err=0.045, m_D_err=0.345)
            #axs[row, col].fill_between(t_values, min_fit_1, max_fit_2, color='orange', alpha=0.5)
            
            axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.0, m_A=1.55, m_D=0.71), color='orange', label='CLAS12 only with fixed D(0)') #*1000*4
            min_fit_1, max_fit_2 = calculate_max_min_fun(cross_section_holo_2, -t_values, Epho_value, D_0=-1.0, m_A=1.55, m_D=0.71, D_0_err=0.101, m_A_err=0.045, m_D_err=0.345)
            axs[row, col].fill_between(t_values, min_fit_1, max_fit_2, color='orange', alpha=0.5)

        if fit_CLAS12_1D:
            print("fit number ",i)
            result_fit = fit_1D(mask_hallB, bin_limits[i], bin_limits[i+1], debug=True)
            
            test_fit_ROOT = fit_1D_Root(mask_hallB, bin_limits[i], bin_limits[i+1], debug=True)
            sigma_0_ROOT = test_fit_ROOT.GetParameter(0)
            m_S_ROOT = test_fit_ROOT.GetParameter(1)
            
            mass_radius = convert_Gev_to_fm * np.sqrt(12)/(result_fit.params['m_S'].value)
            error_mass_radius = (result_fit.params['m_S'].stderr/result_fit.params['m_S'].value)*convert_Gev_to_fm * np.sqrt(12)/(result_fit.params['m_S'].value)
            
            axs[row, col].plot(t_values, model_dipole(-t_values, sigma_0 = result_fit.params['sigma_0'], m_S = result_fit.params['m_S']), label=' LM_fit ')
            axs[row, col].plot(t_values, model_dipole(-t_values, sigma_0 = sigma_0_ROOT, m_S = m_S_ROOT), label=' ROOT fit ')
            axs[row, col].annotate(r'$d\sigma/dt_0$ = ' + f"{result_fit.params['sigma_0'].value:.2f}" + r'$\pm$' + f"{result_fit.params['sigma_0'].stderr:.2f}" + r' nb/GeV$^2$', xy=(1, 0.9), xycoords='axes fraction', ha='right', va='top')
            axs[row, col].annotate(r'$m_S$ = ' + f"{result_fit.params['m_S'].value:.2f}"+ r'$\pm$' + f"{result_fit.params['m_S'].stderr:.2f}" + r' GeV', xy=(1, 0.8), xycoords='axes fraction', ha='right', va='top')
            axs[row, col].annotate(r'$r_m$ = ' + f"{mass_radius:.2f}"+ r'$\pm$' + f"{error_mass_radius:.2f}" + r' fm', xy=(1, 0.7), xycoords='axes fraction', ha='right', va='top')
            
            global Eg_values_CLAS12, eg_errors_up_CLAS12, eg_errors_down_CLAS12, mA_values_CLAS12
            global mA_errors_CLAS12, mass_radius_values_CLAS12, mass_radius_errors_CLAS12
            Eg_values_CLAS12 = np.append(Eg_values_CLAS12, Epho_value)
            eg_errors_up_CLAS12 = np.append(eg_errors_up_CLAS12, bin_limits[i+1]-Epho_value)
            eg_errors_down_CLAS12 = np.append(eg_errors_down_CLAS12, Epho_value-bin_limits[i])
            mA_values_CLAS12 = np.append(mA_values_CLAS12, result_fit.params['m_S'].value) 
            mA_errors_CLAS12 = np.append(mA_errors_CLAS12, result_fit.params['m_S'].stderr) 
            mass_radius_values_CLAS12 = np.append(mass_radius_values_CLAS12, mass_radius) 
            mass_radius_errors_CLAS12 = np.append(mass_radius_errors_CLAS12, error_mass_radius)
            
        # Remove upper and right spines
        axs[row, col].spines['top'].set_visible(False)
        axs[row, col].spines['right'].set_visible(False)
        
       
        if row == (axs.shape[0]-1) and col == (0):
            if(log):
                axs[row, col].legend(loc='lower left', ncol=2, frameon=False, prop={'size': 8})  # Example legend labels
            else:
                axs[row, col].legend(loc='upper right', ncol=2, bbox_to_anchor=(1.0, 0.85), frameon=False, prop={'size': 8})  # Example legend labels

    # Set x-axis label only on the bottom row subplots
    for ax in axs[-1, :]:
        ax.set_xlabel(r'$-t~[GeV^2]$ ')
        
    for i in range(len(bin_limits)-1, len(axs.flatten())):
        row = i // nb_columns
        col = i % nb_columns
        axs[row, col].axis('off')  
        

    fig.text(0.02, 0.5, r'$\frac{d\sigma}{dt}~[nb/GeV^2]$', ha='center', va='center', rotation='vertical')

    plt.tight_layout()
    #plt.show()
    fig.savefig(name)


Plot_Slice()
Plot_Slice(name = "Slice_fit_lin.png", log = False)
Plot_Slice(name="Fit_Clas12.png", nb_columns = 2, bin_limits = [8.2, 9.28, 10.0, 10.6], log = True, only_CLAS12=True, only_Model=False, only_Fit=True)
Plot_Slice(name="Fit_all.png", nb_columns = 4, log = True, only_Model=False, only_Fit=True)
Plot_Slice(name="Data_all.png", nb_columns = 4, log = True, only_Model=False, only_Fit=False)
Plot_Slice(name="Fit_Only_Clas12.png", nb_columns = 2, bin_limits = [8.2, 9.28, 10.0, 10.6], log = False, only_CLAS12=True)

Plot_Slice(name="Data_Only_Clas12.png", nb_columns = 2, bin_limits = [8.2, 9.28, 10.0, 10.6], log = True, only_CLAS12=True, only_Model=False, only_Fit=False, fit_CLAS12_1D=True)


################################################
################   2D plot   ###################
################################################
def TwoDplot():
    fig = plt.figure(figsize=(10, 10))

    plt.scatter(E_gamma_HB, t_HB, color=color_CLAS12, label=label_HallB, s=100)
    plt.scatter(E_gamma_HC, t_HC, color=color_hall_c, label=label_HallC)
    plt.scatter(E_gamma_HD, t_HD, color=color_glueX, label=label_HallD)



    plt.xlabel(r'$E_{\gamma}~[GeV]$')
    plt.ylabel(r'$-t~[GeV^2]$ ')


    # Define the range of W and t values
    Epho_values = np.linspace(7.5, 11.5, 1000)  # Assuming W starts from Mv + Mn = 3.1 + 0.938 GeV
    t_values = np.linspace(-10, 1, 1000)

    # Choose a specific value of xi
    xi_values = [ 0.35, 0.4, 0.5, 0.6, 0.7, 0.8]

    # Plot the curves for each xi
    for xi in xi_values:
        t_points = []
        Epho_points = []
        for t in t_values:
            for Epho in Epho_values:
                if abs(xi_GPD_approach(Epho, t) - xi) < 0.0001 and (t<T_min(Epho) and t>T_max(Epho)) :  # You can adjust the tolerance
                    #print(T_min(Epho))
                    t_points.append(-t)
                    Epho_points.append(Epho)
        if xi == 0.4:          
            plt.plot(Epho_points,t_points, color='b', label=r'Constant $\xi$ line')
        else:
            plt.plot(Epho_points,t_points, color='b')
        plt.annotate(r'$\xi$={}'.format(xi), xy=(Epho_points[0], t_points[0]), xytext=(10, 0), textcoords='offset points')


    Epho_values1 = np.linspace(7.5, 11.5, 1600)
    t_min_values = T_min(Epho_values1)
    t_max_values = T_max(Epho_values1)

    plt.plot(Epho_values1,-t_min_values, color='red', label='Limit of the phase space')
    plt.plot(Epho_values1,-t_max_values, color='red')

    plt.legend(frameon=False, prop={'size': 15})
    #axs[row, col].legend(loc='lower left', ncol=2, frameon=False, prop={'size': 8})  
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().set_ylim(top=10.)
    plt.gca().set_ylim(bottom=0.0)
    plt.gca().set_xlim(right=11.5)
    plt.gca().set_xlim(left=8.2)
    #plt.show()
    fig.savefig("2D.png")

#TwoDplot()

################################################
############### Plot the Form Factors ##########
################################################
fig = plt.figure(figsize=(8, 4))
gs = GridSpec(1, 2, width_ratios=[1, 1])

# Generate t values
min_t = 0.0
max_t = 4.
t_values = np.linspace(0.0, 4.5, 100)

# Compute function values
A_g_values_HB, A_g_upper_HB, A_g_lower_HB= Dipole_g_band(-t_values, 0.414, 0.008, 1.55, 0.04, expo=3)
A_g_values_HC, A_g_upper_HC, A_g_lower_HC= Dipole_g_band(-t_values, 0.414, 0.008, 1.575, 0.059, expo=3)
A_g_values_HC_GPD, A_g_upper_HC_GPD, A_g_lower_HC_GPD= Dipole_g_band(-t_values, 0.414, 0.008, 2.71, 0.19, expo=3)
A_g_values_all, A_g_upper_all, A_g_lower_all= Dipole_g_band(-t_values, 0.414, 0.008, 1.86, 0.057, expo=3)


D_g_values_HB, D_g_upper_HB, D_g_lower_HB = Dipole_g_band(-t_values, -1.8, 0.52, 1.12, 0.21, expo=3)
D_g_values_HC, D_g_upper_HC, D_g_lower_HC = Dipole_g_band(-t_values, -1.8, 0.52, 0.63, 0.23, expo=3)
D_g_values_all, D_g_upper_all, D_g_lower_all = Dipole_g_band(-t_values, -1.01, 0.10, 1.66, 0.18, expo=3)


# Plot the function
ax0 = plt.subplot(gs[0])
ax0.fill_between(t_values, A_g_lower_HC, A_g_upper_HC, color=color_hall_c, alpha=0.5, label='Hall C Holo')
ax0.plot(t_values, A_g_values_HC, color=color_hall_c)

ax0.fill_between(t_values, A_g_lower_HC_GPD, A_g_upper_HC_GPD, color='green', alpha=0.5, label='Hall C GPD old')
ax0.plot(t_values, A_g_values_HC_GPD, color='green')

ax0.fill_between(t_values, A_g_lower_HB, A_g_upper_HB, color=color_CLAS12, alpha=0.5, label='Hall B Holo')
ax0.plot(t_values, A_g_values_HB, color=color_CLAS12)

ax0.fill_between(t_values, A_g_lower_all, A_g_upper_all, color='purple', alpha=0.5, label='Hall B, C, D')
ax0.plot(t_values, A_g_values_all, color='purple')


ax1 = plt.subplot(gs[1])
ax1.fill_between(t_values, D_g_lower_HC, D_g_upper_HC, color=color_hall_c, alpha=0.5, label='Hall C only')
ax1.plot(t_values, D_g_values_HC, color=color_hall_c)

ax1.fill_between(t_values, D_g_lower_HB, D_g_upper_HB, color=color_CLAS12, alpha=0.5, label='Hall C only')
ax1.plot(t_values, D_g_values_HB, color=color_CLAS12)

ax1.fill_between(t_values, D_g_lower_all, D_g_upper_all, color='purple', alpha=0.5, label='Hall C only')
ax1.plot(t_values, D_g_values_all, color='purple')

ax0.set_xlabel('-t')
ax0.set_ylabel(r'$A_g(t)$ ')
ax0.set_ylim(top=0.5)
ax0.set_ylim(bottom=0.0)
ax0.set_xlim(left=min_t)
ax0.set_xlim(right=max_t)


ax1.set_xlabel('-t')
ax1.set_ylabel(r'$D_g(t)$ ')
#ax1.set_yscale('log')
#ax1.set_ylim(bottom=7e-4)
#ax1.set_ylim(top=5)
ax1.set_ylim(top=0.5)
ax1.set_xlim(left=min_t)
ax1.set_xlim(right=max_t)



ax0.legend(frameon=False)
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
fig.savefig("Test_Fit.png")
#plt.show()



### Figure radius###
# Data Hall C
Eg_values_hall_c = np.array([9.175, 9.375, 9.475, 9.625, 9.775, 9.925, 10.075, 10.225, 10.375, 10.525])
eg_errors_hall_c = np.array([0.025, 0.009, 0.017, 0.018, 0.043, 0.038, 0.024, 0.027, 0.033, 0.038])
mA_values_hall_c = np.array([1.69, 2.30, 1.59, 1.52, 1.18, 1.26, 1.43, 1.31, 1.33, 1.32])
mA_errors_hall_c = np.array([0.45, 0.47, 0.22, 0.17, 0.16, 0.17, 0.16, 0.14, 0.19, 0.20])
mass_radius_values_hall_c = np.array([0.404, 0.297, 0.430, 0.450, 0.579, 0.544, 0.477, 0.520, 0.515, 0.517])
mass_radius_errors_hall_c = np.array([0.107, 0.061, 0.059, 0.052, 0.079, 0.075, 0.053, 0.057, 0.074, 0.079])
color_rec_hall_c='#79863c'

# Data CLAS12
#Eg_values_CLAS12 = np.array([9.00,9.86])
#eg_errors_up_CLAS12 = np.array([0.28,0.50])
#eg_errors_down_CLAS12 = np.array([0.80, 0.58])
#mA_values_CLAS12 = np.array([1.43, 1.27])
#mA_errors_CLAS12 = np.array([0.23, 0.09])
#mass_radius_values_CLAS12 = np.array([0.48, 0.54])
#mass_radius_errors_CLAS12 = np.array([0.08, 0.04])

# Data GlueX
Eg_values_GlueX = np.array([8.93, 9.86, 10.82])
eg_errors_up_GlueX = np.array([0.35, 0.5, 0.62])
eg_errors_down_GlueX = np.array([0.73, 0.58, 0.46])
mA_values_GlueX = np.array([1.105, 1.472, 1.313])
mA_errors_GlueX = np.array([0.168, 0.075, 0.049])
mass_radius_values_GlueX = np.array([0.619, 0.464, 0.521])
mass_radius_errors_GlueX = np.array([0.094, 0.024, 0.020])

# Calculate mean and standard deviation for the first four data points of the blue dataset
mean_hall_c_1 = np.mean(mass_radius_values_hall_c[:4])
std_dev_hall_c_1 = np.std(mass_radius_values_hall_c[:4])

mean_hall_c_2 = np.mean(mass_radius_values_hall_c[4:])
std_dev_hall_c_2 = np.std(mass_radius_values_hall_c[4:])

# Create a square grid with two columns for the two plots
fig = plt.figure(figsize=(8, 8))
gs = GridSpec(1, 2, width_ratios=[1, 1])

# Plot mA as a function of Eg with rotated graph
ax0 = plt.subplot(gs[0])
ax0.errorbar(mA_values_hall_c, Eg_values_hall_c, xerr=mA_errors_hall_c,  marker='o', linestyle='', color=color_hall_c, capsize=5, label='Hall C-007 (2022)')
ax0.errorbar(mA_values_GlueX, Eg_values_GlueX, xerr=mA_errors_GlueX, yerr=[eg_errors_down_GlueX, eg_errors_up_GlueX], marker='o', linestyle='', color=color_glueX, capsize=5, label='GlueX (2023)')
ax0.errorbar(mA_values_CLAS12, Eg_values_CLAS12, xerr=mA_errors_CLAS12, yerr=[eg_errors_down_CLAS12, eg_errors_up_CLAS12], marker='o', linestyle='', color=color_CLAS12, capsize=5, label='CLAS12')
ax0.set_ylabel('Eg (GeV)')
ax0.set_xlabel('mA (GeV)')

# Plot mass radius as a function of Eg with rotated graph
ax1 = plt.subplot(gs[1], sharey=ax0)
ax1.errorbar(mass_radius_values_hall_c, Eg_values_hall_c, xerr=mass_radius_errors_hall_c, marker='o', linestyle='', color=color_hall_c, capsize=5, label='Hall C-007')
ax1.errorbar(mass_radius_values_CLAS12, Eg_values_CLAS12, xerr=mass_radius_errors_CLAS12, yerr=[eg_errors_down_CLAS12, eg_errors_up_CLAS12], marker='o', linestyle='', color=color_CLAS12, capsize=5, label='CLAS12')
ax1.errorbar(mass_radius_values_GlueX, Eg_values_GlueX, xerr=mass_radius_errors_GlueX, yerr=[eg_errors_down_GlueX, eg_errors_up_GlueX], marker='o', linestyle='', color=color_glueX, capsize=5, label='GlueX')
ax1.set_xlabel('Mass Radius (fm)')
ax1.tick_params(axis='y', labelleft=False)  # Remove y-axis label on the right plot

# Add rectangle for the mean and standard deviation of the first four data points of the blue dataset
#rect = patches.Rectangle((min(Eg_values[:4]), mean_blue - std_dev_blue / 2), max(Eg_values[:4]) - min(Eg_values[:4]), std_dev_blue, linewidth=1, edgecolor='none', facecolor='blue')
rect_1 = patches.Rectangle((mean_hall_c_1, Eg_values_hall_c[0]),std_dev_hall_c_1, Eg_values_hall_c[3]-Eg_values_hall_c[0], linewidth=1, edgecolor='none', facecolor=color_rec_hall_c, alpha=0.5)
ax1.add_patch(rect_1)

rect_2 = patches.Rectangle((mean_hall_c_2, Eg_values_hall_c[4]),std_dev_hall_c_2, Eg_values_hall_c[9]-Eg_values_hall_c[4], linewidth=1, edgecolor='none', facecolor=color_rec_hall_c, alpha=0.5)
ax1.add_patch(rect_2)

# Add legend to the right plot without box
legend = ax0.legend(loc='upper right', ncol=1, fancybox=True, shadow=False, fontsize=13)
frame = legend.get_frame()
frame.set_linewidth(0)  # Remove the box around the legend

# Remove top and right border lines
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)

# Remove top and right border lines
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Force the y-axis to finish at 10.6
ax0.set_ylim(ax0.get_ylim()[0], 11.5)

# Force the x-axis to finish at 0.7
ax1.set_xlim(ax1.get_xlim()[0], 0.8)

ax0.set_xlabel('m$_s$ (GeV)')  # Replace 'mA' with 'ms'
ax0.set_ylabel('E$γ$ (GeV)')  # Replace 'Eg' with 'Eγ'


plt.text(0.5, 0.5, 'Preliminary', fontsize=100, rotation=45, ha='center', va='center', alpha=0.2, transform=fig.transFigure)

# Adjust layout to prevent overlap and remove space between graphs
plt.subplots_adjust(wspace=0)

# Save the figure as a PNG file
plt.savefig('M_S_and_r_m.png')