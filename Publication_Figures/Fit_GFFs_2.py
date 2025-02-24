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
from scipy.integrate import quad
from scipy.optimize import curve_fit

### DISCLAIMER ###
# This code has become a monster. Use it at your own risk.
##################

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
    
def eta_Holo(Epho,t):
    s = (Mn**2 + 2*Mn*Epho)
    eta = Mv**2 / ( 2*(s-Mn**2) - Mv**2 + t )
    return eta

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

Conv_factor = 0.3894*(10**6)# GeV-2 to nb
    
alpha_S = 0.3
e_q = (2/3)
Phi_V_2 = 1.0952/(4. * math.pi)
print("the factor in front of GPD model ",Conv_factor*( (alpha_EM * e_q**2 /4) * ( (16.0*math.pi*alpha_S)**2 / (3*Mv**3) ) * Phi_V_2 ))

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

def gpd_based_model(input_var,  C_0, m_A, m_C, A_0=0.414):
    t, Epho = input_var
    s = (Mn**2 + 2*Mn*Epho)
    const1_unitconv=70.231
    x = -t
    M=0.938
    mJpsi=3.097
    xi = xi_GPD_approach(Epho, t)#(-x - mJpsi**2)/(2*M**2+mJpsi**2 - 2*s +x)
    A = A_0 / ( 1- (t/(m_A**2)) )**3
    C = C_0 / ( 1- (t/(m_C**2)) )**3
    H2= A + 4* xi**2 *C
    E2= -xi**2 *C*4
    const_2=1/pow(s-M*M,2)
    myfit =const1_unitconv*const_2*(4 /pow(xi, 4))*((1 + x / (4 * M * M)) * E2 *E2 - 2 * E2 * (H2 + E2) + (1 - xi * xi) * (H2 + E2) *(H2+E2))
    return myfit

def gpd_based_model_old(input_var,  C_0, m_A, m_C, A_0=0.414):
    t, Epho = input_var
    s = (Mn**2 + 2*Mn*Epho)
    const1_unitconv=70.231
    x = -t
    M=0.938
    mJpsi=3.097
    xi = (-x - mJpsi**2)/(2*M**2+mJpsi**2 - 2*s +x)
    A = A_0 / ( 1- (t/(m_A**2)) )**3
    C = C_0 / ( 1- (t/(m_C**2)) )**3
    H2= A + 4* xi**2 *C
    E2= -xi**2 *C*4
    const_2=1/pow(s-M*M,2)
    myfit =const1_unitconv*const_2*(1 /pow(xi, 4))*((1 + x / (4 * M * M)) * E2 *E2 - 2 * E2 * (H2 + E2) + (1 - xi * xi) * (H2 + E2) *(H2+E2))
    return myfit



def cross_section_GPD_2(input_var, C_0, m_A, m_C, A_0=0.414):
    
    t, Epho = input_var
    
    M=0.938
    Mj=3.097
    Q=0.0
    const1_unitconv=70.231
    s=M*M-Q*Q+2*M*Epho
    q_gamma = (1/(2*pow(s,0.5)))*pow(s*s-2*(-Q*Q+M*M)*s+(-Q*Q-M*M)*(-Q*Q-M*M),0.5)
    q_jpsi = (1/(2*pow(s,0.5)))*pow(s*s-2*(Mj*Mj+M*M)*s+(Mj*Mj-M*M)*(Mj*Mj-M*M),0.5)
    E_jpsi = pow(Mj*Mj+q_jpsi*q_jpsi,0.5)
    E_gamma = pow(-Q*Q+q_gamma*q_gamma,0.5)
    P_z = (s-M*M)/(2*pow(s,0.5))
    P_0 = pow(P_z*P_z+M*M,0.5)
    P_plus = (P_0+P_z)/(pow(2,0.5))
    P_prime_z = -(-t+Mj*Mj-2*E_gamma*E_jpsi)/(2*q_gamma)
    K=q_jpsi
    P_prime_0=pow((-K)*(-K)+M*M,0.5)
    P_prime_plus=(P_prime_0+P_prime_z)/pow(2,0.5)
    xi = (P_plus-P_prime_plus)/(P_plus+P_prime_plus)
    A_conf=1.0; #0.8; #5/4
    # dipole-dipole
    #H2=A0/pow((1+x/(par[0]*par[0])),2)+pow(2*xi,2)*par[1]/pow((1+x/(par[2]*par[2])),2);
    #E2=-pow(2*xi,2)*par[1]/pow((1+x/(par[2]*par[2])),2);
    #tripole-tripole
    H2=A_0/(A_conf*pow((1-t/(m_A*m_A)),3))+pow(2*xi,2)*C_0/(A_conf*pow((1-t/(m_C*m_C)),3))
    E2=-pow(2*xi,2)*C_0/(A_conf*pow((1-t/(m_C*m_C)),3))
    #dipole-tripole
    #H2=A0/pow((1+x/(par[0]*par[0])),3)+pow(2*xi,2)*par[1]/pow((1+x/(par[2]*par[2])),3);
    #E2=-pow(2*xi,2)*par[1]/pow((1+x/(par[2]*par[2])),3);
    const_2=1/pow(s-M*M,2)
    myfit =const1_unitconv*const_2*(4 /pow(xi, 4))*((1 - t / (4 * M * M)) * E2 *E2 - 2 * E2 * (H2 + E2) + (1 - xi * xi) * (H2 + E2) *(H2+E2))
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

def cross_section_holo_2_C(input_var, C_0, m_A, m_C, A_0=0.414, expo=3, N=2.032): #Formula in PHYSICAL REVIEW D 106, 086004 (2022)
    t, Epho = input_var
    
    N_e = N #in nb GeV-2
    #N_e = 7.768 * 0.389 * 2/3
       
    s = (Mn**2 + 2*Mn*Epho)
    
    tilde_F = ((s-Mn**2)/2)**4
    
    eta = Mv**2 / ( 2*(s-Mn**2) - Mv**2 + t )

    #A_0 = 0.414
    A_t = A_0 / ( 1- (t/(m_A**2)) )**expo
    D_t = 4. * C_0 / ( 1- (t/(m_C**2)) )**expo

    
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

def holographic_model(input_var, m_A, C_0, m_C, A_0=0.414): # This is the model in Nature
    t, Epho = input_var
    A0=0.414
    mJpsi = 3.096
    mp=0.938
    Mn=mp
    s = (Mn**2 + 2*Mn*Epho)
    
    x = -t
    const_h = 1 / (64* np.pi * (s - Mn ** 2) ** 2 * 4  * Mn ** 2)
    c2 = 7.768#5.93#7.0 #6.0 #7.768
    # mS=1
    eta = mJpsi ** 2 / (2 * (s - Mn ** 2) - mJpsi ** 2 - x)
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
        * (2 * x + 8 * Mn ** 2)
    )

def Dipole_g(t, A_0, m_A, expo_D):
    return A_0 / ( 1- (t/(m_A**2)) )**expo_D

def Dipole_g_band(t_values, A_0, A_0_error, m_A, m_A_error, expo=2, n_samples=100):
    """
    Compute the central value, upper, and lower bounds of the Dipole_g function by sampling.
    
    Parameters:
        t_values (array): Input t values.
        A_0 (float): Central value of A_0.
        A_0_error (float): Error range for A_0.
        m_A (float): Central value of m_A.
        m_A_error (float): Error range for m_A.
        expo (int): Exponent for the Dipole_g function.
        n_samples (int): Number of samples per parameter range.
    
    Returns:
        tuple: Central values, upper bounds, lower bounds.
    """
    # Create sampled ranges for A_0 and m_A
    A_0_range = np.linspace(A_0 - A_0_error, A_0 + A_0_error, n_samples)
    m_A_range = np.linspace(m_A - m_A_error, m_A + m_A_error, n_samples)
    
    # Create a grid of parameter combinations
    A_0_grid, m_A_grid = np.meshgrid(A_0_range, m_A_range)
    
    # Flatten the grids for evaluation
    A_0_flat = A_0_grid.ravel()
    m_A_flat = m_A_grid.ravel()
    
    # Evaluate Dipole_g for all combinations
    Dipole_g_values_all = np.array([
        Dipole_g(t_values, a0, ma, expo_D=expo)
        for a0, ma in zip(A_0_flat, m_A_flat)
    ])
    
    # Reshape results back to grid shape for bounds computation
    Dipole_g_values_all = Dipole_g_values_all.reshape(len(A_0_range), len(m_A_range), -1)
    
    # Compute bounds over all sampled combinations
    Dipole_g_upper = np.max(Dipole_g_values_all, axis=(0, 1))
    Dipole_g_lower = np.min(Dipole_g_values_all, axis=(0, 1))
    Dipole_g_central = Dipole_g(t_values, A_0, m_A, expo_D=expo)
    
    return Dipole_g_central, Dipole_g_upper, Dipole_g_lower

#def Dipole_g_band(t_values, A_0, A_0_error, m_A, m_A_error, expo=2):
#    # Compute function values with central parameters
#    Dipole_g_values = Dipole_g(t_values, A_0, m_A, expo_D=expo)
#    
#    # Define all parameter combinations
#    A_0_variations = [A_0 + A_0_error, A_0 - A_0_error]
#    m_A_variations = [m_A + m_A_error, m_A - m_A_error]
#
#    # Compute Dipole_g for all combinations of parameter variations
#    Dipole_g_all_variations = [
#        Dipole_g(t_values, A_0_var, m_A_var, expo_D=expo)
#        for A_0_var in A_0_variations
#        for m_A_var in m_A_variations
#    ]
#    
#    # Compute upper and lower bounds
#    Dipole_g_upper = np.max(Dipole_g_all_variations, axis=0)
#    Dipole_g_lower = np.min(Dipole_g_all_variations, axis=0)
#    
#    return Dipole_g_values, Dipole_g_upper, Dipole_g_lower

def model_expo(t_values, sigma_0, B_G):
    return sigma_0* np.exp(B_G*t_values)

def model_dipole(t_values, sigma_0, m_S):
    return Dipole_g(t_values, sigma_0, m_S, 4)

def integrated_model_GPD( Epho, t1, t2, C_0, m_A, m_C, A_0=0.414):
    # Define the function of t to integrate
    def integrand1(t):
        return gpd_based_model((t, Epho), C_0, m_A, m_C, A_0)

    # Perform numerical integration
    integral, _ = quad(integrand1, t1, t2)
    return integral

def integrated_model_Holo(Epho, t1, t2, C_0, m_A, m_C, A_0=0.414):
    # Define the function of t to integrate
    def integrand2(t):
        return cross_section_holo_2_C((t, Epho), C_0, m_A, m_C, A_0)

    # Perform numerical integration
    integral1, _ = quad(integrand2, t1, t2)
    return integral1

def integrated_model_dipole(t1, t2, s_0, m_S):
    # Define the function of t to integrate
    def integrand2(t):
        return model_dipole(t, s_0, m_S)

    # Perform numerical integration
    integral1, _ = quad(integrand2, t1, t2)
    return integral1

def integrated_model_expo(t1, t2, s_0, B_G):
    # Define the function of t to integrate
    def integrand2(t):
        return model_expo(t, s_0, B_G)

    # Perform numerical integration
    integral1, _ = quad(integrand2, t1, t2)
    return integral1


################################################
################ Integrated CS #################
################################################

def Plot_int_CS():
    plot_all =True
    with open("../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_Latex_Table.txt", "r") as file:
        lines = file.readlines()
        data_HB_int = np.array([[float(value) for value in line.split("&")] for line in lines])

    systematics_CLAS12 = np.array([
        #[130.453, 16.9437, 2.9274, 1.55966, 6.02919, 7.81237, 13.7265, 33.8826, 19.5421, 32.7588],
        #[ 0, 41.6628, 148.573, 54.8454, 34.6713, 134.213, 16.218, 31.8332, 36.9747, 88.9859],
        [26.0444 , 16.9021 , 23.9025 , 9.40902 , 11.8196 , 14.1321 , 11.5646 , 12.6512 , 10.3102 , 20.4096], 
        [26.0444 , 16.9021 , 23.9025 , 9.40902 , 11.8196 , 14.1321 , 11.5646 , 12.6512 , 10.3102 , 20.4096],
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

    
    if plot_all:
        plt.errorbar(data_HD_int[0],data_HD_int[1], yerr=data_HD_int[2], fmt='o', capsize=0, color='#ce0037', label='GlueX (2023)')
    plt.errorbar(data_HB_int[0],data_HB_int[2], xerr=data_HB_int[1], yerr=data_HB_int[3], fmt='o', capsize=0, color='#00b4d8', label='CLAS12 (this work)')
    print(T_max(10), T_min(10))
    Epho_model = np.linspace(8.20, 11.5, 100) 

    integrated_values_Holo = [
        #holographic_model((-t_values, Epho_value), m_A=1.575, C_0=-.45, m_C=1.12),label='Holo Nature')
        #cross_section_holo_2((-t_values, Epho_value), m_A=1.612, D_0=-1.275, m_D=0.9),label='Holo 2022')
        #integrated_model(holographic_model, Epho, T_max(Epho), T_min(Epho),   m_A=1.575, C_0=-.45, m_C=1.12) for Epho in Epho_model
        #integrated_model2(cross_section_holo_2, Epho, T_max(Epho), T_min(Epho),  A_0 = 0.430 , m_A=1.612, D_0=-1.275, m_D=0.9) for Epho in Epho_model
        integrated_model_Holo(Epho, T_max(Epho), T_min(Epho), C_0=-0.32, m_A=1.612, m_C=0.9, A_0 = 0.430) for Epho in Epho_model

    ]

    integrated_values_Holo_up = [integrated_model_Holo(Epho, T_max(Epho), T_min(Epho), C_0=-0.32, m_A=1.612, m_C=0.963, A_0 = 0.430) for Epho in Epho_model]
    integrated_values_Holo_down = [integrated_model_Holo(Epho, T_max(Epho), T_min(Epho), C_0=-0.32, m_A=1.612, m_C=0.963, A_0 = 0.430) for Epho in Epho_model]


    integrated_values_GPD = [
        integrated_model_GPD(Epho, T_max(Epho),  T_min(Epho),  m_A=2.07, C_0=-1.21, m_C=0.91) for Epho in Epho_model #T_max(Epho), T_min(Epho), 
    ]

    integrated_values_GPD_up = [integrated_model_GPD(Epho, T_max(Epho),  T_min(Epho),  m_A=2.07, C_0=-1.21+0.37, m_C=0.91) for Epho in Epho_model]
    integrated_values_GPD_down = [integrated_model_GPD(Epho, T_max(Epho),  T_min(Epho),  m_A=2.07, C_0=-1.21-0.37, m_C=0.91) for Epho in Epho_model]

    #integrated_values_GPD_up = [integrated_model_GPD(Epho, T_max(Epho),  T_min(Epho),  m_A=2.07+0.05, C_0=-1.21+0.37, m_C=0.91-0.1) for Epho in Epho_model]
    #integrated_values_GPD_down = [integrated_model_GPD(Epho, T_max(Epho),  T_min(Epho),  m_A=2.07-0.05, C_0=-1.21-0.37, m_C=0.91+0.1) for Epho in Epho_model]

   
   
    # Plot the results
    if plot_all:
        #plt.plot(Epho_model, integrated_values_GPD, label="G-J-L-Y 2023 (GPD)", color='green')
        #plt.plot(Epho_model, integrated_values_GPD_up, label="G-J-L-Y 2023 (GPD)", color='green')
        #plt.plot(Epho_model, integrated_values_GPD_down, label="G-J-L-Y 2023 (GPD)", color='green')
        plt.fill_between(Epho_model, integrated_values_GPD_down, integrated_values_GPD_up, color='green', alpha=0.5, label="G-J-L-Y 2023 (GPD)")
        plt.plot(Epho_model, integrated_values_Holo, label="M-Z 2022 (Holo.)")

    plt.xlabel(r'$E_{\gamma}~[GeV]$')
    plt.ylabel(r'$\sigma~[nb]$ ')
    
    plt.gca().set_ylim(bottom=0.01)
    plt.gca().set_xlim(left=8.2, right=11.5)

    tot_syst = []

    if plot_all:
        for x, y, yerr, yfill in zip(data_HD_int[0], data_HD_int[1], data_HD_int[2], data_HD_int[3]):
            plt.fill_between([x - 0.03, x + 0.03], y - np.sqrt(yfill**2 + (y*0.2)**2), y + np.sqrt(yfill**2 + (y*0.2)**2), color='#ce0050', alpha=0.3)

    for i, (x, y, syst_low, syst_high) in enumerate(zip(data_HB_int[0], data_HB_int[2], systematics_CLAS12[0], systematics_CLAS12[1])):
        #plt.fill_between([x - 0.03, x + 0.03], y*(1-0.2), y *(1+0.2), color='#00b4d0', alpha=0.3)
        syst_normalization = 16.
        syst_charge = 1.2
        syst_rad_int = [16., 10., 9., 5., 5., 3., 3., 2., 1., 0.]
        error_syst = (np.sqrt(syst_low*syst_low + syst_normalization*syst_normalization + syst_charge*syst_charge + syst_rad_int[i]*syst_rad_int[i])/100.)
        tot_syst.append(error_syst.item())
        plt.fill_between([x - 0.03, x + 0.03], y*(1-error_syst), y *(1+error_syst), color='#00b4d0', alpha=0.3)

    print("Tot. syst. int")
    print(tot_syst)
    print("End of ot. syst. int")
    plt.legend(loc='upper left', ncol=2,  frameon=False, prop={'size': 15})
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    if plot_all:
        plt.gca().set_yscale('log')
        plt.gca().set_ylim(top=10)
    #plt.show()
    fig.savefig("Integrated_CS_only.png", dpi=300)

Plot_int_CS()

def test_integral():

    print("HERE ",cross_section_holo_2_C((-2, 10.45), C_0=-0.32, m_A=1.612, m_C=0.9, A_0 = 0.430))

    Epho = 10.45
    t_values = np.linspace( T_max(Epho), T_min(Epho), 100)
    #integrand_values = [gpd_based_model((t, Epho), m_A=2.07, C_0=-1.21, m_C=0.91) for t in t_values]
    integrand_values = [cross_section_holo_2_C((t, Epho), C_0=-0.32, m_A=1.612, m_C=0.9, A_0 = 0.430) for t in t_values]
    
    plt.figure(figsize=(8, 6))
    plt.plot(t_values, integrand_values, label=f"Integrand for $E_\\gamma = {Epho:.2f}$ GeV")
    plt.xlabel(r"$t$", fontsize=14)
    plt.ylabel("Integrand Value", fontsize=14)
    plt.title("Integrand as a Function of $t$", fontsize=16)
    plt.grid()
    plt.legend(fontsize=12)

    plt.savefig("test_GPD1.png", dpi=300, bbox_inches='tight')  # High resolution, no extra whitespace
#test_integral()

def Plot_int_diff_CS():

    with open("../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_Latex_Table.txt", "r") as file:
        lines = file.readlines()
        data_HB_int = np.array([[float(value) for value in line.split("&")] for line in lines])

    with open("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_8.20_9.28_Latex_Table.txt", "r") as file:
        lines = file.readlines()
        data_HB_1 = np.array([[float(value) for value in line.split("&")] for line in lines])
        bin_1 = [0.77, 1.00, 1.5, 2.0, 2.5, 4.5]
        bin_diff_1 = np.diff(bin_1)
        print(bin_diff_1)
        print(data_HB_1)
        result_1 = np.sum(data_HB_1[4] * bin_diff_1)

    with open("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_9.28_10.00_Latex_Table.txt", "r") as file:
        lines = file.readlines()
        data_HB_2 = np.array([[float(value) for value in line.split("&")] for line in lines])
        bin_2 = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 6.0]
        bin_diff_2 = np.diff(bin_2)
        result_2 = np.sum(data_HB_2[4] * bin_diff_2)

    with open("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_10.00_10.60_Latex_Table.txt", "r") as file:
        lines = file.readlines()
        data_HB_3 = np.array([[float(value) for value in line.split("&")] for line in lines])
       #bin_3 = [0.5, 0.8, 1.1, 1.6, 3.0]
        bin_3 = [0.5, 0.7, 0.9, 1.1, 1.3, 1.6, 2.0, 4.5]
        bin_diff_3 = np.diff(bin_3)
        result_3 = np.sum(data_HB_3[4] * bin_diff_3)

    
    integrated_values_dipole_1_u = integrated_model_dipole(T_max(8.9), T_min(8.9), s_0=1.86-1.70 , m_S=1.37 )
    integrated_values_dipole_1_d = integrated_model_dipole(T_max(8.9), T_min(8.9), s_0=1.86+1.70, m_S=1.37 )
    integrated_values_expo_1_u = integrated_model_expo(T_max(8.9), T_min(8.9), s_0=0.83-0.42 , B_G=1.02 )
    integrated_values_expo_1_d = integrated_model_expo(T_max(8.9), T_min(8.9), s_0=0.83+0.42 , B_G=1.02 )

    integrated_values_dipole_2_u = integrated_model_dipole(T_max(9.65), T_min(9.65), s_0=2.93-1.04 , m_S=1.34)
    integrated_values_dipole_2_d = integrated_model_dipole(T_max(9.65), T_min(9.65), s_0=2.93+1.04 , m_S=1.34)
    integrated_values_expo_2_u = integrated_model_expo(T_max(9.65), T_min(9.65), s_0=1.46-0.28 , B_G=1.08 )
    integrated_values_expo_2_d = integrated_model_expo(T_max(9.65), T_min(9.65), s_0=1.46+0.28 , B_G=1.08 )

    integrated_values_dipole_3_u = integrated_model_dipole(T_max(10.28), T_min(10.28), s_0=3.07-0.77 , m_S=1.23 )
    integrated_values_dipole_3_d = integrated_model_dipole(T_max(10.28), T_min(10.28), s_0=3.07+0.77 , m_S=1.23 )
    integrated_values_expo_3_u = integrated_model_expo(T_max(10.28), T_min(10.28), s_0=1.50-0.29 , B_G=1.31 )
    integrated_values_expo_3_d = integrated_model_expo(T_max(10.28), T_min(10.28), s_0=1.50+0.29 , B_G=1.31 )

    



    fig = plt.figure(figsize=(10, 7))

    
    plt.errorbar(data_HB_int[0],data_HB_int[2], xerr=data_HB_int[1], yerr=data_HB_int[3], fmt='o', capsize=0, color='#00b4d8', label='CLAS12 (this work)')
    #plt.plot([8.2, 9.28], [integrated_values_dipole_1, integrated_values_dipole_1], color='blue', linewidth=2)
    #plt.plot([9.28, 10], [integrated_values_dipole_2, integrated_values_dipole_2], color='blue', linewidth=2)
    #plt.plot([10., 10.6], [integrated_values_dipole_3, integrated_values_dipole_3], yerr=integrated_values_dipole_3_u-integrated_values_dipole_3_d, color='blue', linewidth=2)
   
    #plt.errorbar(10.3, integrated_values_dipole_3, xerr= 0.3, yerr=abs(integrated_values_dipole_3_u-integrated_values_dipole_3_d), color='blue', linewidth=2)
    
    plt.fill_between([8.2, 9.28], integrated_values_expo_1_u, integrated_values_expo_1_d, color='blue', alpha=0.3, label='Integrated exponential model')
    plt.fill_between([9.28, 10.0], integrated_values_expo_2_u, integrated_values_expo_2_d, color='blue', alpha=0.3)
    plt.fill_between([10.0,10.6], integrated_values_expo_3_u, integrated_values_expo_3_d, color='blue', alpha=0.3)

    plt.fill_between([8.2, 9.28], integrated_values_dipole_1_u, integrated_values_dipole_1_d, color='orange', alpha=0.3, label='Integrated dipole model')
    plt.fill_between([9.28, 10.0], integrated_values_dipole_2_u, integrated_values_dipole_2_d, color='orange', alpha=0.3)
    plt.fill_between([10.0,10.6], integrated_values_dipole_3_u, integrated_values_dipole_3_d, color='orange', alpha=0.3)

   
   
    plt.xlabel(r'$E_{\gamma}~[GeV]$')
    plt.ylabel(r'$\sigma~[nb]$ ')
    
    plt.gca().set_ylim(bottom=0.01)
    plt.gca().set_xlim(right=11.5)

    
    plt.legend(loc='upper left', ncol=2,  frameon=False, prop={'size': 15})
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    plt.gca().set_yscale('log')
    plt.gca().set_ylim(top=10)
    #plt.show()
    fig.savefig("Integrated_diff_CS.png")

Plot_int_diff_CS()


################################################
############## Differential  CS ################
################################################
def Plot_diff_CS(input, input_rad_syst, array, name):
    with open(input, "r") as file:
        lines = file.readlines()
        data_HB_int = np.array([[float(value) for value in line.split("&")] for line in lines])

    systematics_CLAS12 = np.array([
        array, 
        array,
    ])

    fig = plt.figure(figsize=(10, 7))

    plt.errorbar(data_HB_int[0],data_HB_int[4], xerr=data_HB_int[1], yerr=data_HB_int[5], fmt='o', capsize=0, color='#00b4d8')#, label='CLAS12 (this work)')

    plt.xlabel(r'$-t~[GeV^2]$ ')
    plt.ylabel(r'$\frac{d\sigma}{dt}~[nb/GeV^2]$')
    plt.gca().set_ylim(top=2)
    plt.gca().set_ylim(bottom=0.01)
    plt.gca().set_xlim(left=0.0)
    plt.gca().set_xlim(right=5.0)

    tot_syst = []

    for x, y, syst_low, syst_high in zip(data_HB_int[0], data_HB_int[4], systematics_CLAS12[0], systematics_CLAS12[1]):
        #plt.fill_between([x - 0.03, x + 0.03], y*(1-0.2), y *(1+0.2), color='#00b4d0', alpha=0.3)
        syst_normalization = 16.
        syst_charge = 1.2
        syst_rad = input_rad_syst
        print(x,y)
        error_syst = (np.sqrt(syst_low*syst_low + syst_normalization*syst_normalization + syst_charge*syst_charge +syst_rad*syst_rad)/100.)
        tot_syst.append(float(error_syst.item()))
        plt.fill_between([x - 0.1, x + 0.1], y*(1-error_syst), y *(1+error_syst), color='#00b4d0', alpha=0.3)

    print("Tot. syst. ", input)
    print(tot_syst)
    plt.legend(loc='upper right', ncol=2,  frameon=False, prop={'size': 15})
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().set_yscale('log')
    #plt.show()
    fig.savefig(name+".png")

Plot_diff_CS("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_8.20_9.28_Latex_Table.txt", 8, [32.7321 , 8.42031 , 18.5855 , 26.4685 , 24.9405], "bin1")
Plot_diff_CS("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_9.28_10.00_Latex_Table.txt", 4, [25.534 , 23.4012 , 8.19891 , 19.1962 , 9.90248 , 15.459 , 36.8322 , 29.3957 , 9.09458], "bin2")
Plot_diff_CS("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_10.00_10.60_Latex_Table.txt", 1, [33.5969 , 25.5195 , 10.8375 , 13.3547 , 14.3888 , 32.07 , 20.3599], "bin3")

################################################
################ CLAS12 ########################
################################################
with open("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_8.20_9.28_Latex_Table.txt", "r") as file:
    lines = file.readlines()
    data_1 = np.array([[float(value) for value in line.split("&")] for line in lines])
     
with open("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_9.28_10.00_Latex_Table.txt", "r") as file:
    lines = file.readlines()
    data_2 = np.array([[float(value) for value in line.split("&")] for line in lines]) 
    
with open("../CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_10.00_10.60_Latex_Table_test.txt", "r") as file:
    lines = file.readlines()
    data_3 = np.array([[float(value) for value in line.split("&")] for line in lines]) 

# Concatenate data from both files
#t_HB = np.concatenate((data_1[0], data_2[0]))#, data_3[0]))
#E_gamma_HB = np.concatenate((data_1[2], data_2[2]))#, data_3[2]))
#dsdt_nb_HB = np.concatenate((data_1[4], data_2[4]))#, data_3[4]))
#err_dsdt_nb_HB = np.concatenate((data_1[5], data_2[5]))#, data_3[5]))

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

#I know this is the worse posible way to store data, I will do it in df in another life


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

t_data = np.concatenate((t_data_C, t_data_B, t_data_D))
E_data = np.concatenate((E_data_C, E_data_B, E_data_D))
cs_data = np.concatenate((cs_data_C, cs_data_B, cs_data_D))
sigma_cs = np.concatenate((sigma_cs_C, sigma_cs_B, sigma_cs_D))
xi_all = [float(xi_GPD_approach(Epho, t)) for Epho, t in zip(E_data, t_data)]

#t_data = np.concatenate(( t_data_B, t_data_D))
#E_data = np.concatenate(( E_data_B, E_data_D))
#cs_data = np.concatenate(( cs_data_B, cs_data_D))
#sigma_cs = np.concatenate(( sigma_cs_B, sigma_cs_D))
#xi_all = [float(xi_GPD_approach(Epho, t)) for Epho, t in zip(E_data, t_data)]

print(xi_all)
cut_xi=0.40
mask_xi = [xi > cut_xi for xi in xi_all]

t_data_GPD = t_data[mask_xi]
E_data_GPD = E_data[mask_xi]
cs_data_GPD = cs_data[mask_xi]
sigma_cs_GPD = sigma_cs[mask_xi]


#t_data_C_B = np.concatenate((t_data_D, t_data_B))
#E_data_C_B = np.concatenate((E_data_D, E_data_B))
#cs_data_C_B = np.concatenate((cs_data_D, cs_data_B))
#sigma_cs_C_B = np.concatenate((sigma_cs_D, sigma_cs_B))

#t_data_C_B = np.concatenate((t_data_C, t_data_B))
#E_data_C_B = np.concatenate((E_data_C, E_data_B))
#cs_data_C_B = np.concatenate((cs_data_C, cs_data_B))
#sigma_cs_C_B = np.concatenate((sigma_cs_C, sigma_cs_B))


#t_data_C_B = np.concatenate((t_data_C, t_data_D))
#E_data_C_B = np.concatenate((E_data_C, E_data_D))
#cs_data_C_B = np.concatenate((cs_data_C, cs_data_D))
#sigma_cs_C_B = np.concatenate((sigma_cs_C, sigma_cs_D))

t_data_C_B = t_data_B
E_data_C_B = E_data_B
cs_data_C_B = cs_data_B
sigma_cs_C_B = sigma_cs_B

xi_C_B = [float(xi_GPD_approach(Epho, t)) for Epho, t in zip(E_data_C_B, t_data_C_B)]

print(xi_all)
mask_xi_C_B = [xi > cut_xi for xi in xi_C_B]
#mask_bug = [ (E>4.0 and E<9.55) or (E>9.7 and E<10.15) and E>10.30 for E, t in zip(E_data_C_B, t_data_C_B)]
#mask_xi_C_B = np.logical_and(mask_xi_C_B, mask_bug)

t_data_C_B_GPD = t_data_C_B[mask_xi_C_B]
E_data_C_B_GPD = E_data_C_B[mask_xi_C_B]
cs_data_C_B_GPD = cs_data_C_B[mask_xi_C_B]
sigma_cs_C_B_GPD = sigma_cs_C_B[mask_xi_C_B]

#mask_C_B_holo = [E>4.0 or (E>10.15 and E<10.30 and t>-1.0) for E,t in zip(E_data_C_B, t_data_C_B)]
#print(t_data_C_B)
#t_data_C_B_Holo = t_data_C_B[mask_C_B_holo]
#E_data_C_B_Holo = E_data_C_B[mask_C_B_holo]
#cs_data_C_B_Holo = cs_data_C_B[mask_C_B_holo]
#sigma_cs_C_B_Holo = sigma_cs_C_B[mask_C_B_holo]
#
#print((E_data_C_B_Holo, t_data_C_B_Holo))


experiment_array = np.concatenate((np.full_like(E_data_C, 1), np.full_like(E_data_B, 2),np.full_like(E_data_D, 3)))
#experiment_array = np.full_like(E_data, 2)


model_GPD = Model(gpd_based_model)
params_GPD = model_GPD.make_params(C_0 = -.45, m_A = 1.5, m_C = 1.0, debug=False)
params_GPD['m_C'].min = 0
#params_GPD['C_0'].min = -1.0
params_GPD['m_A'].min = 0
params_GPD['A_0'].vary = False
#params_GPD['debug'].vary = False

params_GPD_B = model_GPD.make_params(C_0 = -1.478, m_A = 2.07, m_C = 0.91, debug=False)
params_GPD_B['m_C'].min = 0
#params_GPD['C_0'].min = -1.0
params_GPD_B['m_A'].min = 0
params_GPD_B['A_0'].vary = False
params_GPD_B['m_C'].vary = False
#params_GPD_B['C_0'].vary = False
fix_C0=False
fix_mA=False
fix_mC=True
#params_GPD['debug'].vary = False


#model_holo = Model(cross_section_holo_2)
#params_holo = model_holo.make_params(D_0=-1.8, m_A=1.612, m_D=0.963)
#params_holo['m_A'].min = 0
#params_holo['D_0'].min = -3
##params_holo['D_0'].vary = False
#params_holo['A_0'].vary = False
#params_holo['N'].vary = False
#params_holo['expo'].vary = False

model_holo = Model(cross_section_holo_2_C)
params_holo = model_holo.make_params(m_A=1.5, C_0=-.45, m_C=1.25)
#params_holo['m_A'].min = 0
#params_holo['C_0'].min = -3
#params_holo['D_0'].vary = False
params_holo['A_0'].vary = False
params_holo['N'].vary = False
params_holo['expo'].vary = False

params_holo_B = model_holo.make_params(m_A=1.641, C_0=-.45, m_C=1.12)
#params_holo['m_A'].min = 0
#params_holo['C_0'].min = -3
#params_holo['D_0'].vary = False
params_holo_B['A_0'].vary = False
params_holo_B['m_C'].vary = False
params_holo_B['N'].vary = False
params_holo_B['expo'].vary = False


apply_xi_cut = True
if apply_xi_cut:
    print("Do the GPD fit C and B")
    result_GPD_C_B = model_GPD.fit(cs_data_C_B_GPD, input_var = (t_data_C_B_GPD, E_data_C_B_GPD), params=params_GPD_B, weights=1/sigma_cs_C_B_GPD)
    print(result_GPD_C_B.fit_report())
    print(result_GPD_C_B.params['m_C'].value)
    #print(result_GPD_C_B.covar[0][2]/np.sqrt(result_GPD_C_B.covar[0][0]*result_GPD_C_B.covar[2][2]))
    print("///////////////////\n")

    print("///////////////////")
    print("Results all data GPD")
    result_GPD = model_GPD.fit(cs_data_GPD, input_var = (t_data_GPD, E_data_GPD), params=params_GPD, weights=1/sigma_cs_GPD)
    print(result_GPD.fit_report())
    print(result_GPD.covar[0][2]/np.sqrt(result_GPD.covar[0][0]*result_GPD.covar[2][2]))
    print("///////////////////\n")

    if False:
        # Fit the model and cross check with scipy
        popt, pcov = curve_fit(
            gpd_based_model,
            (t_data_C_B_GPD.ravel(), E_data_C_B_GPD.ravel()),
            cs_data_C_B_GPD.ravel(),
            p0=[-.45, 1.5, 1.25],
            sigma=sigma_cs_C_B_GPD.ravel(),
            absolute_sigma=True,
        )

        # Extract fitted parameters and uncertainties
        fitted_params = popt
        uncertainties = np.sqrt(np.diag(pcov))

        print("Fitted parameters for GPD -  C and B:")
        param_names = ["C0", "mA", "mC"]
        for name, value, uncertainty in zip(param_names, fitted_params, uncertainties):
            print(f"{name}: {value:.3f} ± {uncertainty:.3f}")

else:
    print("///////////////////")
    print("Do the GPD fit C and B")
    result_GPD_C_B = model_GPD.fit(cs_data_C_B, input_var = (t_data_C_B, E_data_C_B), params=params_GPD, weights=1/sigma_cs_C_B)
    print(result_GPD_C_B.fit_report())
    print(result_GPD_C_B.params['m_C'].value)
    print("///////////////////\n")

    print("///////////////////")
    print("Results all data GPD")
    result_GPD = model_GPD.fit(cs_data, input_var = (t_data, E_data), params=params_GPD, weights=1/sigma_cs)
    print(result_GPD.fit_report())
    print("///////////////////\n")

#print("///////////////////")
#print("Do the GPD fit C and B")
#result_GPD_C_B = model_GPD.fit(cs_data_C_B, input_var = (t_data_C_B, E_data_C_B), params=params_GPD, weights=1/sigma_cs_C_B)
#print(result_GPD_C_B.fit_report())
#print(result_GPD_C_B.params['m_C'].value)
#print("///////////////////\n")
#
#print("///////////////////")
#print("Results all data GPD")
#result_GPD = model_GPD.fit(cs_data, input_var = (t_data, E_data), params=params_GPD, weights=1/sigma_cs)
#print(result_GPD.fit_report())
#print("///////////////////\n")




print("///////////////////")
print("Do the Holo fit C and B")
#mask = (E_data_C_B>9.7) | (E_data_C_B<9.55)
#result_Holo_C_B = model_holo.fit(cs_data_C_B[mask], input_var = (t_data_C_B[mask], E_data_C_B[mask]), params=params_holo, weights=1/sigma_cs_C_B[mask])
result_Holo_C_B = model_holo.fit(cs_data_C_B, input_var = (t_data_C_B, E_data_C_B), params=params_holo_B, weights=1/sigma_cs_C_B)#, method="lbfgsb")
print(result_Holo_C_B.fit_report())
print(result_Holo_C_B.params['m_C'].value)
#print(result_Holo_C_B.covar[0][2]/np.sqrt(result_Holo_C_B.covar[0][0]*result_Holo_C_B.covar[2][2]))
print("///////////////////\n")

# Fit the model and cross check with scipy
popt, pcov = curve_fit(
    cross_section_holo_2_C,
    (t_data_C_B.ravel(), E_data_C_B.ravel()),
    cs_data_C_B.ravel(),
    p0=[-.45, 1.5, 1.25],
    sigma=sigma_cs_C_B.ravel(),
    absolute_sigma=True,
)

# Extract fitted parameters and uncertainties
fitted_params = popt
uncertainties = np.sqrt(np.diag(pcov))

print("Fitted parameters:")
param_names = ["C0", "mA", "mC"]
for name, value, uncertainty in zip(param_names, fitted_params, uncertainties):
    print(f"{name}: {value:.3f} ± {uncertainty:.3f}")


print("///////////////////")
print("Results all data holo")
result_holo = model_holo.fit(cs_data, input_var = (t_data, E_data), params=params_holo, weights=1/sigma_cs)
print(result_holo.fit_report())
print(result_holo.covar[0][2]/np.sqrt(result_holo.covar[0][0]*result_holo.covar[2][2]))
print(result_holo.covar[0][1]/np.sqrt(result_holo.covar[0][0]*result_holo.covar[1][1]))
print("///////////////////\n")



#print("///////////////////")
#print("Results C")
#result_GPD_C = model_GPD.fit(cs_data_C, input_var = (t_data_C, E_data_C), params=params_GPD, weights=1/sigma_cs_C)
#print(result_GPD_C.fit_report())
#print("///////////////////\n")
#
#print("///////////////////")
#print("Results C and D")
#result_holo_D_C = model_holo.fit(cs_data_D_C, input_var = (t_data_D_C, E_data_D_C), params=params_holo, weights=1/sigma_cs_D_C)
#print(result_holo_D_C.fit_report())
#print("///////////////////\n")

#print("Results B")
#params_holo['D_0'].value = -1.01331202
#params_holo['D_0'].vary = False
#result_holo_B = model_holo.fit(cs_data_B, input_var = (t_data_B, E_data_B), params=params_holo, weights=1/sigma_cs_B)
#print(result_holo_B.fit_report())
#print("///////////////////\n")

#print("Results D")
#params_holo['D_0'].vary = True
#result_holo_D = model_holo.fit(cs_data_D, input_var = (t_data_D, E_data_D), params=params_holo, weights=1/sigma_cs_D)
#print(result_holo_D.fit_report())
#print("///////////////////\n")



################################################
### Do the 1D fits #############################
################################################
model_expo_fit = Model(model_expo)
params_expo = model_expo_fit.make_params(sigma_0=2.0, B_G=1.87)
def fit_expo(mask_hallB, E_1, E_2, debug=False):
    print(params_expo)
    mask = (E_data >= E_1) & (E_data <= E_2)
    result_expo = model_expo_fit.fit(cs_data[mask][mask_hallB], t_values=t_data[mask][mask_hallB], params=params_expo, weights=1/sigma_cs[mask][mask_hallB])
    
    if debug:
        print(result_expo.fit_report())
        
    return result_expo

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
#########   Define the 1D fit table here   #####
################################################
W_values_CLAS12 = np.array([]) 
BG_CLAS12 = np.array([])
BG_errors_CLAS12 = np.array([])

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
def three_d_fit():
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
def Plot_Slice(name = "Slice_fit.png", nb_columns = 3, bin_limits = [8.20, 9.10, 9.25, 9.40, 9.55, 9.70, 9.85, 10.0, 10.15, 10.30, 10.45, 10.6, 11.5], log = True, only_CLAS12=False, only_Model = False, only_Fit = False, fit_CLAS12_1D = False, fit_expo_1D = False):#, t_data=t_data, E_data=E_data, cs_data=cs_data, sigma_cs=sigma_cs):
    t_data = np.concatenate((t_data_C, t_data_B, t_data_D))
    E_data = np.concatenate((E_data_C, E_data_B, E_data_D))
    cs_data = np.concatenate((cs_data_C, cs_data_B, cs_data_D))
    sigma_cs = np.concatenate((sigma_cs_C, sigma_cs_B, sigma_cs_D))

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

        tmin=-T_min(bin_limits[i]*1.1)-0.1
        
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
        axs[row, col].annotate(r'$E_{\gamma}$ in ' + f"[{bin_limits[i]:.2f}, {bin_limits[i+1]:.2f}] GeV", xy=(1, 1), xycoords='axes fraction', ha='right', va='top', fontsize=14)
        axs[row, col].set_ylabel(' ')

        t_values = np.linspace(tmin,  8.0+0.2, 100)
        if(bin_limits[i]!=8.2):  
            t_values = np.linspace(-T_min(bin_limits[i]),  8.0+0.2, 100)
        Epho_value = np.mean(E_display) #(bin_limits[i] + (bin_limits[i+1]-bin_limits[i])/2)
        s_values = np.full_like(t_values, Mn**2 + 2*Mn*Epho_value)
        

        #axs[row, col].plot(t_values, cross_section_GPD((-t_values, Epho_value), C_0 = -0.48, m_A=1.64, m_C=1.07),label='GPD N=0.1') #1000*4
        
        if only_Model:

            if(bin_limits[i]!=8.2):
                axs[row, col].axvline(x=-T_min(bin_limits[i]), color='g', linestyle='--')
            #axs[row, col].plot(t_values, cross_section_holo((-t_values, Epho_value), m_A=1.5),label='Holo') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.04, m_A=1.87, m_D=1.67),label='Holo All') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.8, m_A=1.575, m_D=1.12),label='Holo Nature') #*1000*4

            axs[row, col].plot(t_values, holographic_model((-t_values, Epho_value), m_A=1.575, C_0=-.45, m_C=1.12),label='Nature 2023 (Holo.)') #*1000*4
            #axs[row, col].plot(t_values, holographic_model((-t_values, Epho_value), m_A=1.59, C_0=-.38, m_C=1.206),label='Holo All') #*1000*4
            axs[row, col].plot(t_values, gpd_based_model_old((-t_values, Epho_value), m_A=2.71, C_0=-0.2, m_C=1.28),label='Nature 2023 (GPD) ') #*1000*4
            axs[row, col].plot(t_values, gpd_based_model((-t_values, Epho_value), m_A=2.07, C_0=-1.21, m_C=0.91),label='G-J-L-Y 2023 (GPD)') #*1000*4


            #axs[row, col].plot(t_values, gpd_based_model((-t_values, Epho_value), m_A=result_GPD.params['m_A'].value, C_0=result_GPD.params['C_0'].value, m_C=result_GPD.params['m_C'].value),label='GPD all') #*1000*4
            #axs[row, col].plot(t_values, gpd_based_model((-t_values, Epho_value), m_A=result_GPD_C.params['m_A'].value, C_0=result_GPD_C.params['C_0'].value, m_C=result_GPD_C.params['m_C'].value),label='GPD C new') #*1000*4

            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), A_0 = 0.430, m_A=1.612, D_0=-1.275, m_D=0.9),label='Holo 2022')
            axs[row, col].plot(t_values, cross_section_holo_2_C((-t_values, Epho_value),  A_0 = 0.430 , m_A=1.612, C_0=-0.32, m_C=0.9),label='M-Z 2022 (Holo.)') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), m_A=1.59, D_0=-.38*4, m_D=1.206),label='Holo 2022 with Nature param') #*1000*4

            print(result_GPD_C_B.params['m_A'].value, result_GPD_C_B.params['C_0'].value, result_GPD_C_B.params['m_C'].value)

            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.93, m_A=1.641, m_D=1.7, A_0=0.429, expo=3, N=2.549),label='Tripole Lattice') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-10, m_A=1.13, m_D=0.48, A_0=0.58, expo=2, N=3.244),label='Dipole Lattice') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.8, m_A=1.575, m_D=1.12),label='HC Holo fit 2') #*1000*4
            #axs[row, col].plot(t_values, holographic_model((-t_values, s_values), m_A=1.575, C_0=-.45, m_C=1.12),label='HC Holo fit') #*1000*4 
            #axs[row, col].plot(t_values, cross_section_GPD_1((-t_values, Epho_value), C_0=-0.45, m_A=1.575, m_C=1.12, A_0=0.414),label='HC GPD - Holo fit') #*1000*4
            #axs[row, col].plot(t_values, cross_section_GPD_1((-t_values, Epho_value), C_0=-0.2, m_A=1.7, m_C=1.28, A_0=0.414),label='HC GPD - GPD fit ') #*1000*4

            #axs[row, col].plot(t_values, cross_section_GPD_1((-t_values, Epho_value), A_0=0.414, C_0=-1.763, m_A=1.84, m_C=0.74),label='GPD fit ') #*1000*4
            #axs[row, col].plot(t_values, cross_section_GPD_1((-t_values, Epho_value), A_0=0.414, C_0=-0.45, m_A=1.575, m_C=1.12),label='GPD fit 1 ') #*1000*4
        
        if only_Fit:
            if(bin_limits[i]!=8.2):
                axs[row, col].axvline(x=-T_min(bin_limits[i]), color='g', linestyle='--')
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=result_holo.params['D_0'], m_A=result_holo.params['m_A'], m_D=result_holo.params['m_D']),label='Fit Hall B') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.19, m_A=1.87, m_D=1.56),label=' My HC Holo fit') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.8, m_A=1.55, m_D=0.63),label=' My HB Holo fit with fixed D') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-3.0, m_A=1.55, m_D=0.57),label=' My HB Holo fit with fixed D') #*1000*4
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-0.7, m_A=1.78, m_D=1.54),label=' My HB Holo fit without bin 3') #*1000*4

            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.05, m_A=2.04, m_D=1.96),label=' My HD Holo fit') #*1000*4

            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.0, m_A=1.86, m_D=1.66), color='green', label='Combined fit') #*1000*4
            #min_fit, max_fit = calculate_max_min_fun(cross_section_holo_2, -t_values, Epho_value, D_0=-1.0, m_A=1.86, m_D=1.66, D_0_err=0.101, m_A_err=0.057, m_D_err=0.176)
            #axs[row, col].fill_between(t_values, min_fit, max_fit, color='green', alpha=0.5)
            
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.0, m_A=1.55, m_D=0.71),label='CLAS12 only with fixed D(0)') #*1000*4
            #min_fit_1, max_fit_2 = calculate_max_min_fun(cross_section_holo_2, -t_values, Epho_value, D_0=-1.0, m_A=1.55, m_D=0.71, D_0_err=0.101, m_A_err=0.045, m_D_err=0.345)
            #axs[row, col].fill_between(t_values, min_fit_1, max_fit_2, color='orange', alpha=0.5)
            
            #axs[row, col].plot(t_values, cross_section_holo_2((-t_values, Epho_value), D_0=-1.0, m_A=1.57, m_D=0.84), color='orange', label='CLAS12 only with fixed D(0)') #*1000*4
            #min_fit_1, max_fit_2 = calculate_max_min_fun(cross_section_holo_2, -t_values, Epho_value, D_0=-1.0, m_A=1.57, m_D=0.84, D_0_err=0.101, m_A_err=0.06, m_D_err=0.33)
            #axs[row, col].fill_between(t_values, min_fit_1, max_fit_2, color='orange', alpha=0.5)

            axs[row, col].plot(t_values, cross_section_holo_2_C((-t_values, Epho_value), C_0=result_Holo_C_B.params['C_0'].value, m_A=result_Holo_C_B.params['m_A'].value, m_C=result_Holo_C_B.params['m_C'].value), color='orange', label=rf"Fit Holo. - B only - $\chi^2_\nu = {result_Holo_C_B.redchi:.3f}$")

            axs[row, col].plot(t_values, gpd_based_model((-t_values, Epho_value), C_0=result_GPD_C_B.params['C_0'].value, m_A=result_GPD_C_B.params['m_A'].value, m_C=result_GPD_C_B.params['m_C'].value), color='blue', label=rf"Fit GPD - B only - $\chi^2_\nu = {result_GPD_C_B.redchi:.3f}$")

            axs[row, col].plot(t_values, cross_section_holo_2_C((-t_values, Epho_value), C_0=result_holo.params['C_0'].value, m_A=result_holo.params['m_A'].value, m_C=result_holo.params['m_C'].value), color='red', label=rf"Fit Holo. - all - $\chi^2_\nu = {result_holo.redchi:.3f}$")
            
            axs[row, col].plot(t_values, gpd_based_model((-t_values, Epho_value), C_0=result_GPD.params['C_0'].value, m_A=result_GPD.params['m_A'].value, m_C=result_GPD.params['m_C'].value), color='green', label=rf"Fit GPD - all - $\chi^2_\nu = {result_GPD.redchi:.3f}$")

            

        if fit_CLAS12_1D:
            print("fit number ",i)
            result_fit = fit_1D(mask_hallB, bin_limits[i], bin_limits[i+1], debug=True)
            
            test_fit_ROOT = fit_1D_Root(mask_hallB, bin_limits[i], bin_limits[i+1], debug=True)
            sigma_0_ROOT = test_fit_ROOT.GetParameter(0)
            m_S_ROOT = test_fit_ROOT.GetParameter(1)
            
            mass_radius = convert_Gev_to_fm * np.sqrt(12)/(result_fit.params['m_S'].value)
            error_mass_radius = (result_fit.params['m_S'].stderr/result_fit.params['m_S'].value)*convert_Gev_to_fm * np.sqrt(12)/(result_fit.params['m_S'].value)
            
            #axs[row, col].plot(t_values, model_dipole(-t_values, sigma_0 = result_fit.params['sigma_0'], m_S = result_fit.params['m_S']), label=' LM_fit ')
            axs[row, col].plot(t_values, model_dipole(-t_values, sigma_0 = sigma_0_ROOT, m_S = m_S_ROOT), label=' Dipole fit ')
            axs[row, col].annotate(r'$d\sigma/dt_0$ = ' + f"{result_fit.params['sigma_0'].value:.2f}" + r'$\pm$' + f"{result_fit.params['sigma_0'].stderr:.2f}" + r' nb/GeV$^2$', xy=(1, 0.9), xycoords='axes fraction', ha='right', va='top')
            axs[row, col].annotate(r'$m_S$ = ' + f"{result_fit.params['m_S'].value:.2f}"+ r'$\pm$' + f"{result_fit.params['m_S'].stderr:.2f}" + r' GeV', xy=(1, 0.8), xycoords='axes fraction', ha='right', va='top')
            axs[row, col].annotate(r'$r_m$ = ' + f"{mass_radius:.2f}"+ r'$\pm$' + f"{error_mass_radius:.2f}" + r' fm', xy=(1, 0.7), xycoords='axes fraction', ha='right', va='top')
            axs[row, col].annotate(r'$\chi_\nu^2$ = ' + f"{result_fit.redchi:.2f}", xy=(1, 0.6), xycoords='axes fraction', ha='right', va='top')

            global Eg_values_CLAS12, eg_errors_up_CLAS12, eg_errors_down_CLAS12, mA_values_CLAS12
            global mA_errors_CLAS12, mass_radius_values_CLAS12, mass_radius_errors_CLAS12
            Eg_values_CLAS12 = np.append(Eg_values_CLAS12, Epho_value)
            eg_errors_up_CLAS12 = np.append(eg_errors_up_CLAS12, bin_limits[i+1]-Epho_value)
            eg_errors_down_CLAS12 = np.append(eg_errors_down_CLAS12, Epho_value-bin_limits[i])
            mA_values_CLAS12 = np.append(mA_values_CLAS12, result_fit.params['m_S'].value) 
            mA_errors_CLAS12 = np.append(mA_errors_CLAS12, result_fit.params['m_S'].stderr) 
            mass_radius_values_CLAS12 = np.append(mass_radius_values_CLAS12, mass_radius) 
            mass_radius_errors_CLAS12 = np.append(mass_radius_errors_CLAS12, error_mass_radius)

        if fit_expo_1D:
            print("fit number ",i)
            result_fit = fit_expo(mask_hallB, bin_limits[i], bin_limits[i+1], debug=True)
            
            B_G = result_fit.params['B_G'].value
            sigma_0 = result_fit.params['sigma_0'].value
            
            #axs[row, col].plot(t_values, model_dipole(-t_values, sigma_0 = result_fit.params['sigma_0'], m_S = result_fit.params['m_S']), label=' LM_fit ')
            axs[row, col].plot(t_values, model_expo(-t_values, sigma_0 = result_fit.params['sigma_0'].value, B_G = result_fit.params['B_G'].value), label=' Exponential fit ')
            axs[row, col].annotate(r'$d\sigma/dt_0$ = ' + f"{result_fit.params['sigma_0'].value:.2f}" + r'$\pm$' + f"{result_fit.params['sigma_0'].stderr:.2f}" + r' nb', xy=(1, 0.9), xycoords='axes fraction', ha='right', va='top')
            axs[row, col].annotate(r'$B_G$ = ' + f"{result_fit.params['B_G'].value:.2f}"+ r'$\pm$' + f"{result_fit.params['B_G'].stderr:.2f}" + r' GeV$^{-2}$', xy=(1, 0.8), xycoords='axes fraction', ha='right', va='top')
            axs[row, col].annotate(r'$\chi_\nu^2$ = ' + f"{result_fit.redchi:.2f}", xy=(1, 0.7), xycoords='axes fraction', ha='right', va='top')

            
            global W_values_CLAS12, BG_CLAS12, BG_errors_CLAS12
            W = np.sqrt(0.938**2 + 2*0.938*Epho_value)
            W_values_CLAS12 = np.append(W_values_CLAS12, W)
            BG_CLAS12 = np.append(BG_CLAS12, result_fit.params['B_G'].value) 
            BG_errors_CLAS12 = np.append(BG_errors_CLAS12, result_fit.params['B_G'].stderr) 
            
        # Remove upper and right spines
        axs[row, col].spines['top'].set_visible(False)
        axs[row, col].spines['right'].set_visible(False)
        
       
        if row == (axs.shape[0]-1) and col == (0):
            if(log):
                axs[row, col].legend(loc='lower left', ncol=1, frameon=False, prop={'size': 8})  # Example legend labels
            else:
                axs[row, col].legend(loc='upper right', ncol=1, bbox_to_anchor=(1.0, 0.85), frameon=False, prop={'size': 8})  # Example legend labels

    # Set x-axis label only on the bottom row subplots
    for ax in axs[-1, :]:
        ax.set_xlabel(r'$-t~[GeV^2]$ ', fontsize=14)
        
    for i in range(len(bin_limits)-1, len(axs.flatten())):
        row = i // nb_columns
        col = i % nb_columns
        axs[row, col].axis('off')  
        

    fig.text(0.02, 0.5, r'$\frac{d\sigma}{dt}~[nb/GeV^2]$', ha='center', va='center', rotation='vertical')#, fontsize=14)

    plt.tight_layout()
    #plt.show()
    fig.savefig(name)


#Plot_Slice()
#Plot_Slice(name = "Slice_fit_lin.png", log = False)
Plot_Slice(name="Fit_Clas12.png", nb_columns = 2, bin_limits = [8.2, 9.28, 10.0, 10.6], log = True, only_CLAS12=True, only_Model=False, only_Fit=True)
Plot_Slice(name="Fit_all.png", nb_columns = 4, log = True, only_Model=False, only_Fit=True)
Plot_Slice(name="Data_all.png", nb_columns = 4, log = True, only_Model=False, only_Fit=False)
Plot_Slice(name="Fit_Only_Clas12.png", nb_columns = 2, bin_limits = [8.2, 9.28, 10.0, 10.6], log = False, only_CLAS12=True)

Plot_Slice(name="Fit_GPD.png", nb_columns = 4, only_Model=True, log = True)
Plot_Slice(name="Data_Only_Clas12.png", nb_columns = 2, bin_limits = [8.2, 9.28, 10.0, 10.6], log = True, only_CLAS12=True, only_Model=False, only_Fit=False, fit_CLAS12_1D=True)
Plot_Slice(name="Fit_expo_Only_Clas12.png", nb_columns = 2, bin_limits = [8.2, 9.28, 10.0, 10.6], log = True, only_CLAS12=True, only_Model=False, only_Fit=False, fit_expo_1D=True)



################################################
################   2D plot   ###################
################################################
def TwoDplot():
    fig = plt.figure(figsize=(10, 10))

    plt.scatter(E_gamma_HB, t_HB, color=color_CLAS12, label=label_HallB)#, s=100)
    plt.scatter(E_gamma_HC, t_HC, color=color_hall_c, label=label_HallC)
    plt.scatter(E_gamma_HD, t_HD, color=color_glueX, label=label_HallD)



    plt.xlabel(r'$E_{\gamma}~[GeV]$')
    plt.ylabel(r'$-t~[GeV^2]$ ')


    # Define the range of W and t values
    Epho_values = np.linspace(7.5, 11.5, 1000)  # Assuming W starts from Mv + Mn = 3.1 + 0.938 GeV
    t_values = np.linspace(-10, 1, 1000)

    # Choose a specific value of xi
    xi_values = [ 0.35, 0.4, 0.5, 0.6, 0.7, 0.8]#[ 0.3, 0.35, 0.38, 0.40, 0.45]#

    # Plot the curves for each xi
    for xi in xi_values:
        t_points = []
        Epho_points = []
        for t in t_values:
            for Epho in Epho_values:
                #print(eta_Holo(Epho,t))
                if abs(xi_GPD_approach(Epho, t) - xi) < 0.0001 and (t<T_min(Epho) and t>T_max(Epho)) :  # You can adjust the tolerance
                #if abs(eta_Holo(Epho,t) - xi) < 0.0001 and (t<T_min(Epho) and t>T_max(Epho)) :  # You can adjust the tolerance
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

TwoDplot()

################################################
############### Plot the Form Factors ##########
################################################

def compute_uncertainty_bands(
    f, x0, y0, sigma_x, sigma_y, rho, z_vals, num_points=500, h=1e-5
):
    """
    Compute the uncertainty bands for a function using numerical differentiation.

    Parameters:
        f (callable): The function of one variable z and two parameters (x, y).
        x0 (float): Central value of parameter x.
        y0 (float): Central value of parameter y.
        sigma_x (float): Uncertainty in x.
        sigma_y (float): Uncertainty in y.
        rho (float): Correlation coefficient between x and y.
        z_range (tuple): Range of z values as (z_min, z_max).
        num_points (int): Number of points to calculate in the z range.
        h (float): Step size for numerical differentiation.

    Returns:
        tuple: Arrays of z values, lower uncertainty band, and upper uncertainty band.
    """
    def numerical_derivative(f, var_idx, z, x, y, h=1e-10):
        """Numerically approximate the partial derivative of f."""
        if var_idx == 0:  # Derivative with respect to x
            return (f(z, x + h, y) - f(z, x - h, y)) / (2 * h)
        elif var_idx == 1:  # Derivative with respect to y
            return (f(z, x, y + h) - f(z, x, y - h)) / (2 * h)
        else:
            raise ValueError("Invalid var_idx. Use 0 for x and 1 for y.")

    # Compute f(z) at central values
    f_central = np.array([f(z, x0, y0) for z in z_vals])

    # Propagate uncertainties
    sigma_f = np.sqrt(
        np.array([
            (numerical_derivative(f, 0, z, x0, y0, h) * sigma_x) ** 2
            for z in z_vals
        ]) +
        np.array([
            (numerical_derivative(f, 1, z, x0, y0, h) * sigma_y) ** 2
            for z in z_vals
        ]) +
        2 * rho * np.array([
            numerical_derivative(f, 0, z, x0, y0, h) *
            numerical_derivative(f, 1, z, x0, y0, h) *
            sigma_x * sigma_y
            for z in z_vals
        ])
    )

    # Compute upper and lower bounds
    f_upper = f_central + sigma_f
    f_lower = f_central - sigma_f

    return f_lower, f_upper

def model_tripole(t_values, sigma_0, m_S):
    return Dipole_g(t_values, sigma_0, m_S, 3)

def Plot_GFFs(GPD=False):
    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(1, 2, width_ratios=[1, 1])

    # Generate t values
    min_t = 0.0
    max_t = 4.5
    t_values = np.linspace(0.0, 4.5, 100)

    # Compute function values
    
    #A_g_values_HC, A_g_upper_HC, A_g_lower_HC= Dipole_g_band(-t_values, A_0=0.414, A_0_error= 0.008, m_A=1.575, m_A_error=0.059, expo=3)
    #D_g_values_HC, D_g_upper_HC, D_g_lower_HC = Dipole_g_band(-t_values, A_0=-4*0.45, A_0_error= 4*0.132, m_A=1.12, m_A_error=0.21, expo=3)

    #A_g_values_HC_GPD, A_g_upper_HC_GPD, A_g_lower_HC_GPD= Dipole_g_band(-t_values, A_0=0.414, A_0_error= 0.008, m_A=2.71, m_A_error=0.19, expo=3)
    #D_g_values_HC_GPD, D_g_upper_HC_GPD, D_g_lower_HC_GPD = Dipole_g_band(-t_values, A_0=-4.*0.20, A_0_error= 4*0.11, m_A=1.28, m_A_error=0.5, expo=3)
#
    #A_g_values_C_B_GPD, A_g_upper_C_B_GPD, A_g_lower_C_B_GPD = Dipole_g_band(-t_values, A_0=result_GPD_C_B.params['A_0'].value, A_0_error= 0.008, m_A=result_GPD_C_B.params['m_A'].value, m_A_error=result_GPD_C_B.params['m_A'].stderr, expo=3)
    #D_g_values_C_B_GPD, D_g_upper_C_B_GPD, D_g_lower_C_B_GPD = Dipole_g_band(-t_values, A_0=4.*result_GPD_C_B.params['C_0'].value, A_0_error= 4.*result_GPD_C_B.params['C_0'].stderr, m_A=result_GPD_C_B.params['m_C'].value, m_A_error=result_GPD_C_B.params['m_C'].stderr, expo=3)
#
    #A_g_values_C_B_Holo, A_g_upper_C_B_Holo, A_g_lower_C_B_Holo= Dipole_g_band(-t_values, A_0=result_Holo_C_B.params['A_0'].value, A_0_error= 0.008, m_A=result_Holo_C_B.params['m_A'].value, m_A_error=result_Holo_C_B.params['m_A'].stderr, expo=3)
    #D_g_values_C_B_Holo, D_g_upper_C_B_Holo, D_g_lower_C_B_Holo = Dipole_g_band(-t_values, A_0=4.*result_Holo_C_B.params['C_0'].value, A_0_error= 4.*result_Holo_C_B.params['C_0'].stderr, m_A=result_Holo_C_B.params['m_C'].value, m_A_error=result_Holo_C_B.params['m_C'].stderr, expo=3)
#
#
    #A_g_values_all_GPD, A_g_upper_all_GPD, A_g_lower_all_GPD = Dipole_g_band(-t_values, A_0=result_GPD.params['A_0'].value, A_0_error= 0.008, m_A=result_GPD.params['m_A'].value, m_A_error=result_GPD.params['m_A'].stderr, expo=3)
    #D_g_values_all_GPD, D_g_upper_all_GPD, D_g_lower_all_GPD = Dipole_g_band(-t_values, A_0=4.*result_GPD.params['C_0'].value, A_0_error= 4.*result_GPD.params['C_0'].stderr, m_A=result_GPD.params['m_C'].value, m_A_error=result_GPD.params['m_C'].stderr, expo=3)
#
    #A_g_values_all_Holo, A_g_upper_all_Holo, A_g_lower_all_Holo= Dipole_g_band(-t_values, A_0=result_holo.params['A_0'].value, A_0_error= 0.008, m_A=result_holo.params['m_A'].value, m_A_error=result_holo.params['m_A'].stderr, expo=3)
    #D_g_values_C_B_Holo, D_g_upper_all_Holo, D_g_lower_all_Holo = Dipole_g_band(-t_values, A_0=4.*result_holo.params['C_0'].value, A_0_error= 4.*result_holo.params['C_0'].stderr, m_A=result_holo.params['m_C'].value, m_A_error=result_holo.params['m_C'].stderr, expo=3)

    A_g_upper_HC, A_g_lower_HC = compute_uncertainty_bands(model_tripole, x0 = 0.414, y0 = 1.575, sigma_x=0.008, sigma_y=0.059, rho=0, z_vals=-t_values)
    D_g_upper_HC, D_g_lower_HC = compute_uncertainty_bands(model_tripole, x0 = -4*0.45, y0 = 1.12, sigma_x=4*0.132, sigma_y=0.21, rho=0.9, z_vals=-t_values, h=1e-8)

    print("hahaha", result_GPD_C_B.params['m_A'])

    
    if(fix_mC):
        corr_GPD_B = 0.0
        corr_Holo_B = 0.0
    #if(not fix_C0 and not fix_mA):
    #    corr=result_GPD_C_B.covar[0][2]/np.sqrt(result_GPD_C_B.covar[0][0]*result_GPD_C_B.covar[2][2])
    #elif (fix_mA):
    #    corr=result_GPD_C_B.covar[0][1]/np.sqrt(result_GPD_C_B.covar[1][1]*result_GPD_C_B.covar[0][0])

    A_g_upper_C_B_GPD, A_g_lower_C_B_GPD = compute_uncertainty_bands(model_tripole, x0 = 0.414, y0 =result_GPD_C_B.params['m_A'].value, sigma_x=0.008, sigma_y=result_GPD_C_B.params['m_A'].stderr, rho=0, z_vals=-t_values)
    D_g_upper_C_B_GPD, D_g_lower_C_B_GPD = compute_uncertainty_bands(model_tripole, x0 = 4.*result_GPD_C_B.params['C_0'].value, y0 =0.91, sigma_x=4.*result_GPD_C_B.params['C_0'].stderr, sigma_y=0.10, rho=corr_GPD_B, z_vals=-t_values)

    A_g_upper_C_B_Holo, A_g_lower_C_B_Holo= compute_uncertainty_bands(model_tripole, x0 = 0.414, y0 =result_Holo_C_B.params['m_A'].value, sigma_x=0.008, sigma_y=result_Holo_C_B.params['m_A'].stderr, rho=0, z_vals=-t_values)
    D_g_upper_C_B_Holo, D_g_lower_C_B_Holo = compute_uncertainty_bands(model_tripole, x0 = 4.*result_Holo_C_B.params['C_0'].value, y0 =1.12, sigma_x=4.*result_Holo_C_B.params['C_0'].stderr, sigma_y=0.21, rho=corr_Holo_B, z_vals=-t_values)

    A_g_upper_all_GPD, A_g_lower_all_GPD = compute_uncertainty_bands(model_tripole, x0 = 0.414, y0 =result_GPD.params['m_A'].value, sigma_x=0.008, sigma_y=result_GPD.params['m_A'].stderr, rho=0, z_vals=-t_values)
    D_g_upper_all_GPD, D_g_lower_all_GPD = compute_uncertainty_bands(model_tripole, x0 = 4.*result_GPD.params['C_0'].value, y0 =result_GPD.params['m_C'].value, sigma_x=4.*result_GPD.params['C_0'].stderr, sigma_y=result_GPD.params['m_C'].stderr, rho=result_GPD.covar[0][2]/np.sqrt(result_GPD.covar[0][0]*result_GPD.covar[2][2]), z_vals=-t_values)

    A_g_upper_all_Holo, A_g_lower_all_Holo = compute_uncertainty_bands(model_tripole, x0 = 0.414, y0 =result_holo.params['m_A'].value, sigma_x=0.008, sigma_y=result_holo.params['m_A'].stderr, rho=0, z_vals=-t_values)
    D_g_upper_all_Holo, D_g_lower_all_Holo = compute_uncertainty_bands(model_tripole, x0 = 4.*result_holo.params['C_0'].value, y0 =result_holo.params['m_C'].value, sigma_x=4.*result_holo.params['C_0'].stderr, sigma_y=result_holo.params['m_C'].stderr, rho=result_holo.covar[0][2]/np.sqrt(result_holo.covar[0][0]*result_holo.covar[2][2]), z_vals=-t_values)

    A_g_upper_HallC_GPD, A_g_lower_HallC_GPD = compute_uncertainty_bands(model_tripole, x0 = 0.414, y0 =1.825, sigma_x=0.008, sigma_y=0.063, rho=0, z_vals=-t_values)
    D_g_upper_HallC_GPD, D_g_lower_HallC_GPD = compute_uncertainty_bands(model_tripole, x0 = 4.*-1.478, y0 =0.774, sigma_x=4.*0.467, sigma_y=0.092, rho=0.99, z_vals=-t_values)


    # Plot the function
    if(GPD):
        ax0 = plt.subplot(gs[0])
        ax0.plot([], [], ' ', label=r'$\mathbf{GPD\ Model}$') 
        #ax0.fill_between(t_values, A_g_lower_HC_GPD, A_g_upper_HC_GPD, color=color_hall_c, alpha=0.5, label='Hall C only (Nature)')
        ax0.fill_between(t_values, A_g_lower_C_B_GPD, A_g_upper_C_B_GPD, color='purple', alpha=0.5, label=r'Hall B ($A(0)$ and $m_C$ fixed)')
        ax0.fill_between(t_values, A_g_lower_HallC_GPD, A_g_upper_HallC_GPD, color=color_hall_c, alpha=0.5, label=r'Hall C ($A(0)$ fixed)')
        ax0.fill_between(t_values, A_g_lower_all_GPD, A_g_upper_all_GPD, color='blue', alpha=0.5, label='All data combined ($A(0)$ fixed)')

        ax1 = plt.subplot(gs[1])
        #ax1.fill_between(t_values, D_g_lower_HC_GPD, D_g_upper_HC_GPD, color=color_hall_c, alpha=0.5, label='Hall C only')
        ax1.fill_between(t_values, D_g_lower_C_B_GPD, D_g_upper_C_B_GPD, color='purple', alpha=0.5, label=r'Hall B ($A(0)$ and $m_C$ fixed)')
        ax1.fill_between(t_values, D_g_lower_HallC_GPD, D_g_upper_HallC_GPD, color=color_hall_c, alpha=0.5, label=r'Hall C ($A(0)$ fixed)')
        ax1.fill_between(t_values, D_g_lower_all_GPD, D_g_upper_all_GPD, color='blue', alpha=0.5, label=r'All data combined ($A(0)$ fixed)')

    else:
        ax0 = plt.subplot(gs[0])
        ax0.plot([], [], ' ', label=r'$\mathbf{Holographic\ Model}$') 
        ax0.fill_between(t_values, A_g_lower_C_B_Holo, A_g_upper_C_B_Holo, color='purple', alpha=0.5, label=r'Hall B ($A(0)$ and $m_C$ fixed)')
        ax0.fill_between(t_values, A_g_lower_HC, A_g_upper_HC, color=color_hall_c, alpha=0.5, label=r'Hall C ($A(0)$ fixed)')
        ax0.fill_between(t_values, A_g_lower_all_Holo, A_g_upper_all_Holo, color='blue', alpha=0.5, label=r'All data combined ($A(0)$ fixed)')

        ax1 = plt.subplot(gs[1])
        ax1.fill_between(t_values, D_g_lower_C_B_Holo, D_g_upper_C_B_Holo, color='purple', alpha=0.5, label=r'Hall B ($A(0)$ and $m_C$ fixed)')
        ax1.fill_between(t_values, D_g_lower_HC, D_g_upper_HC, color=color_hall_c, alpha=0.5, label=r'Hall C ($A(0)$ fixed)')
        ax1.fill_between(t_values, D_g_lower_all_Holo, D_g_upper_all_Holo, color='blue', alpha=0.5, label=r'All data combined ($A(0)$ fixed)')


    ax0.set_xlabel(r'-t [GeV$^2$]')
    ax0.set_ylabel(r'$A_g(t)$ ')
    ax0.set_ylim(top=0.5)
    ax0.set_ylim(bottom=0.01)
    ax0.set_xlim(left=min_t)
    ax0.set_xlim(right=max_t)
    #ax0.set_yscale('log')


    ax1.set_xlabel(r'-t [GeV$^2$]')
    ax1.set_ylabel(r'$D_g(t)$ ')
    ax1.set_yscale('symlog', linthresh=0.001)
    #ax1.set_yscale('log')
    #ax1.set_ylim(bottom=7e-4)
    #ax1.set_ylim(top=5)
    ax1.set_ylim(top=-0.001)
    ax1.set_xlim(left=min_t)
    ax1.set_xlim(right=max_t)



    ax0.legend(frameon=False)
    #ax1.legend(frameon=False)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    if(GPD):
        fig.savefig("GFF_GPD.png", dpi=500)
    else:
        fig.savefig("GFF_Holo.png", dpi=500)
    #plt.show()

Plot_GFFs(GPD=False)
Plot_GFFs(GPD=True)



def compute_mass_radius_and_uncertainty(
    m_A, sigma_m_A, A_0, sigma_A_0, C_0, sigma_C_0,
    rho_mA_A0, rho_mA_C0, rho_A0_C0
):
    # Constants
    proton_mass = 0.938

    # Function for the radius
    def radius_func(m_A, A_0, C_0):
        return np.sqrt((18.0 / (m_A**2)) - (6.0 / A_0) * (C_0 / (proton_mass**2)))

    # Partial derivatives
    def dR_dmA(m_A, A_0, C_0):
        return -18.0 / (m_A**3 * radius_func(m_A, A_0, C_0))

    def dR_dA0(m_A, A_0, C_0):
        return (3.0 * C_0) / (A_0**2 * proton_mass**2 * radius_func(m_A, A_0, C_0))

    def dR_dC0(m_A, A_0, C_0):
        return -3.0 / (A_0 * proton_mass**2 * radius_func(m_A, A_0, C_0))

    # Compute partial derivatives at central values
    dR_dmA_val = dR_dmA(m_A, A_0, C_0)
    dR_dA0_val = dR_dA0(m_A, A_0, C_0)
    dR_dC0_val = dR_dC0(m_A, A_0, C_0)

    # Uncertainty propagation formula
    sigma_R = np.sqrt(
        (dR_dmA_val * sigma_m_A)**2 +
        (dR_dA0_val * sigma_A_0)**2 +
        (dR_dC0_val * sigma_C_0)**2 +
        2 * rho_mA_A0 * dR_dmA_val * dR_dA0_val * sigma_m_A * sigma_A_0 +
        2 * rho_mA_C0 * dR_dmA_val * dR_dC0_val * sigma_m_A * sigma_C_0 +
        2 * rho_A0_C0 * dR_dA0_val * dR_dC0_val * sigma_A_0 * sigma_C_0
    )

    # Compute the radius
    radius_val = radius_func(m_A, A_0, C_0)

    return radius_val*convert_Gev_to_fm, sigma_R*convert_Gev_to_fm


def compute_scalar_radius_and_uncertainty(
    m_A, sigma_m_A, A_0, sigma_A_0, C_0, sigma_C_0,
    rho_mA_A0, rho_mA_C0, rho_A0_C0
):
    # Constants
    proton_mass = 0.938

    # Function for the radius
    def radius_func(m_A, A_0, C_0):
        return np.sqrt(
            (18.0 / (m_A**2)) - (18.0 / A_0) * (C_0 / proton_mass**2)
        )

    # Partial derivatives
    def dR_dmA(m_A, A_0, C_0):
        return -18.0 / (m_A**3 * radius_func(m_A, A_0, C_0))

    def dR_dA0(m_A, A_0, C_0):
        return (9.0 * C_0) / (A_0**2 * proton_mass**2 * radius_func(m_A, A_0, C_0))

    def dR_dC0(m_A, A_0, C_0):
        return -9.0 / (A_0 * proton_mass**2 * radius_func(m_A, A_0, C_0))

    # Compute partial derivatives at central values
    dR_dmA_val = dR_dmA(m_A, A_0, C_0)
    dR_dA0_val = dR_dA0(m_A, A_0, C_0)
    dR_dC0_val = dR_dC0(m_A, A_0, C_0)

    # Uncertainty propagation formula
    sigma_R = np.sqrt(
        (dR_dmA_val * sigma_m_A)**2 +
        (dR_dA0_val * sigma_A_0)**2 +
        (dR_dC0_val * sigma_C_0)**2 +
        2 * rho_mA_A0 * dR_dmA_val * dR_dA0_val * sigma_m_A * sigma_A_0 +
        2 * rho_mA_C0 * dR_dmA_val * dR_dC0_val * sigma_m_A * sigma_C_0 +
        2 * rho_A0_C0 * dR_dA0_val * dR_dC0_val * sigma_A_0 * sigma_C_0
    )

    # Compute the radius
    radius_val = radius_func(m_A, A_0, C_0)

    return radius_val*convert_Gev_to_fm, sigma_R*convert_Gev_to_fm

# Example data for 6 models
models = [
    #GPD
    {"name": "All data \n ($A(0)$ fixed)", "A_0": result_GPD.params['A_0'].value, "A_0_error": 0.008, 
     "m_A": result_GPD.params['m_A'].value, "m_A_error": result_GPD.params['m_A'].stderr, 
     "C_0": result_GPD.params['C_0'].value, "C_0_error": result_GPD.params['C_0'].stderr, 
     "m_C": result_GPD.params['m_C'].value, "m_C_error": result_GPD.params['m_C'].stderr,
     "rho_mA_A0": 0.0,
     "rho_mA_C0": result_GPD.covar[0][1]/np.sqrt(result_GPD.covar[0][1]*result_GPD.covar[0][1]),
     "rho_A0_C0": 0.0},

    #{"name": "Hall B and D \n combined", "A_0": result_GPD_C_B.params['A_0'].value, "A_0_error": 0.008, 
    # "m_A": result_GPD_C_B.params['m_A'].value, "m_A_error": result_GPD_C_B.params['m_A'].stderr, 
    # "C_0": result_GPD_C_B.params['C_0'].value, "C_0_error": result_GPD_C_B.params['C_0'].stderr, 
    # "m_C": result_GPD_C_B.params['m_C'].value, "m_C_error": result_GPD_C_B.params['m_C'].stderr,
    # "rho_mA_A0": 0.0,
    # "rho_mA_C0": result_GPD_C_B.covar[0][1]/np.sqrt(result_GPD_C_B.covar[0][1]*result_GPD_C_B.covar[0][1]),
    # "rho_A0_C0": 0.0
    # },

     {"name": "Hall B\n($A(0)$ and \n$m_C$ fixed)", "A_0": 0.414, "A_0_error": 0.008, 
    "m_A": result_GPD_C_B.params['m_A'].value, "m_A_error": result_GPD_C_B.params['m_A'].stderr, 
     "C_0": result_GPD_C_B.params['C_0'].value, "C_0_error": result_GPD_C_B.params['C_0'].stderr, 
     "m_C": 0.91, "m_C_error": 0.10,
     "rho_mA_A0": 0.0,
     "rho_mA_C0": result_GPD_C_B.covar[0][1]/np.sqrt(result_GPD_C_B.covar[0][1]*result_GPD_C_B.covar[0][1]),
     "rho_A0_C0": 0.0},

    {"name": "Hall C\n($A(0)$ fixed)", "A_0": 0.414, "A_0_error": 0.008, 
     "m_A": 1.825, "m_A_error": 0.063, 
     "C_0": -1.478, "C_0_error": 0.467, 
     "m_C": 0.774, "m_C_error": 0.092,
     "rho_mA_A0": 0.0,
     "rho_mA_C0": 0.9,
     "rho_A0_C0": 0.0},

    {"name": "GPD model\n($A(0)$ fixed)", "A_0": 0.414, "A_0_error": 0.008, 
     "m_A": 2.07, "m_A_error": 0.05, 
     "C_0": -1.21, "C_0_error": 0.37, 
     "m_C": 0.91, "m_C_error": 0.10,
     "rho_mA_A0": 0.0,
     "rho_mA_C0": 0.0,
     "rho_A0_C0": 0.0},


    #Holo
    {"name": "Model 6", "A_0": result_holo.params['A_0'].value, "A_0_error": 0.008, 
     "m_A": result_holo.params['m_A'].value, "m_A_error": result_holo.params['m_A'].stderr, 
     "C_0": result_holo.params['C_0'].value, "C_0_error": result_holo.params['C_0'].stderr, 
     "m_C": result_holo.params['m_C'].value, "m_C_error": result_holo.params['m_C'].stderr,
     "rho_mA_A0": 0.0,
     "rho_mA_C0": result_holo.covar[0][1]/np.sqrt(result_holo.covar[0][1]*result_holo.covar[0][1]),
     "rho_A0_C0": 0.0},

    {"name": "Model 5", "A_0": result_Holo_C_B.params['A_0'].value, "A_0_error": 0.008, 
     "m_A": result_Holo_C_B.params['m_A'].value, "m_A_error": result_Holo_C_B.params['m_A'].stderr, 
     "C_0": result_Holo_C_B.params['C_0'].value, "C_0_error": result_Holo_C_B.params['C_0'].stderr, 
     "m_C": 1.12, "m_C_error": 0.21,
     "rho_mA_A0": 0.0,
     "rho_mA_C0": result_Holo_C_B.covar[0][1]/np.sqrt(result_Holo_C_B.covar[0][1]*result_Holo_C_B.covar[0][1]),
     "rho_A0_C0": 0.0},

    {"name": "Model 4", "A_0": 0.414, "A_0_error": 0.008, 
     "m_A": 1.575, "m_A_error": 0.059, 
     "C_0": -0.45, "C_0_error": 0.132, 
     "m_C": 1.12, "m_C_error": 0.21,
     "rho_mA_A0": 0.0,
     "rho_mA_C0": 0.9,
     "rho_A0_C0": 0.0},

    {"name": "GPD model", "A_0": -100000., "A_0_error": 0, 
     "m_A": 0.001, "m_A_error": 0.0, 
     "C_0": -1.20, "C_0_error": 0.0, 
     "m_C": 0.91, "m_C_error": 0.5,
     "rho_mA_A0": 0.0,
     "rho_mA_C0": 0.0,
     "rho_A0_C0": 0.0},
]

# Compute derivatives and errors at t = 0 and t = 1
results = [
    {
        "name": model["name"],
        "derivatives": [
            compute_mass_radius_and_uncertainty(model["m_A"], model["m_A_error"], model["A_0"], model["A_0_error"], model["C_0"], model["C_0_error"], model["rho_mA_A0"], model["rho_mA_C0"], model["rho_A0_C0"]),
            compute_scalar_radius_and_uncertainty(model["m_A"], model["m_A_error"], model["A_0"], model["A_0_error"], model["C_0"], model["C_0_error"], model["rho_mA_A0"], model["rho_mA_C0"], model["rho_A0_C0"]),
        ]
    }
    for model in models
]

print(results)

# Separate into two groups for plotting
left_models = results[:4]
right_models = results[4:]

# Plotting
fig, axs = plt.subplots(1, 2, figsize=(10, 6), sharey=True, gridspec_kw={"wspace": 0})

# Adjust spacing between points
offset = 0.2

# Plot left panel
for i, model in enumerate(left_models):
    for j, (derivative, error) in enumerate(model["derivatives"]):
        y_pos = 2 * i + (-offset if j == 0 else offset)
        axs[0].errorbar(
            derivative, y_pos, xerr=error, fmt='o',
            color='blue' if j == 0 else '#00004D', capsize=5, label=model["name"] if j == 0 else None
        )
        print("Radius", y_pos, derivative, error)
        if (i > 2):
            axs[0].text(
            derivative, y_pos + 0.1, 'Mass' if j == 0 else 'Scalar',
            fontsize=8, ha="center", color='blue' if j == 0 else '#00004D'
            )
axs[0].set_title("GPD model")
axs[0].set_xlabel("Radius [fm]")
axs[0].spines['top'].set_visible(False)
axs[0].spines['right'].set_visible(False)

# Plot right panel
for i, model in enumerate(right_models):
    for j, (derivative, error) in enumerate(model["derivatives"]):
        y_pos = 2 * i + (-offset if j == 0 else offset)
        if model["name"] == "Model 4":
            if j==0:
                derivative=0.755
                error=0.067
            if j==1:
                derivative=1.069
                error=0.126
        axs[1].errorbar(
            derivative, y_pos, xerr=error, fmt='o',
            color='orange' if j == 0 else '#CC6600', capsize=5, label=model["name"] if j == 0 else None
        )
        print("Radius", y_pos, derivative, error)
        if (i > 1):
            axs[1].text(
            derivative, y_pos + 0.1, 'Mass' if j == 0 else 'Scalar',
            fontsize=8, ha="center", color='orange' if j == 0 else '#CC6600'
            )
axs[1].set_title("Holographic model")
axs[1].set_xlabel("Radius [fm]")
axs[1].spines['top'].set_visible(False)
axs[1].spines['right'].set_visible(False)

# Shared y-axis adjustments
axs[0].spines['left'].set_visible(True)
axs[1].spines['left'].set_visible(False)
axs[0].set_yticks([2 * i for i in range(len(left_models))])
axs[0].set_yticklabels([m["name"] for m in left_models])
axs[1].tick_params(axis='y', left=False)
axs[0].set_xlim(left=0)
axs[0].set_xlim(right=3.0)
axs[1].set_xlim(left=0)
axs[1].set_xlim(right=1.6)

labels = axs[0].get_xticklabels()

# Remove the last label
labels[-1].set_visible(False)

# Add central spine between panels
for ax in axs:
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)

# Adjust layout
#fig.tight_layout()
fig.savefig("Plot_radii.png", dpi=300)


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
fig = plt.figure(figsize=(12, 8))
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
ax1.set_xlabel('Mass Radius (fm)', fontsize=16)
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
ax1.set_xlim(ax1.get_xlim()[0], 0.9)
ax1.axvline(x=0.84, color='blue', linewidth=2)
ax1.text(0.84, 9.5, 'CODATA charge radius of the proton', color='blue', fontsize=12, rotation=90, va='center', ha='right')


ax0.set_xlabel('m$_s$ (GeV)', fontsize=16)  # Replace 'mA' with 'ms'
ax0.set_ylabel('E$γ$ (GeV)', fontsize=16)  # Replace 'Eg' with 'Eγ'


plt.text(0.5, 0.5, 'Preliminary', fontsize=100, rotation=45, ha='center', va='center', alpha=0.2, transform=fig.transFigure)

# Adjust layout to prevent overlap and remove space between graphs
plt.subplots_adjust(wspace=0)

# Save the figure as a PNG file
plt.savefig('M_S_and_r_m.png', dpi=300)



def w_to_tau(W) :
    M_jpsi = 3.097
    mp=0.938
    tau = M_jpsi**2/(W**2-mp**2)
    return tau





########## FIGURE BG #######
### Figure radius###
# Data H1 2005
W_H1_2005 = np.array([   44.8, 54.8, 64.8, 74.8, 84.9, 94.9, 104.9, 119.5, 144.1, 180.6, 250.7])
BG_H1_2005 = np.array([  4.13, 4.30, 4.57, 4.46, 4.45, 4.72, 4.79, 4.71, 4.95, 5.08, 5.41])
e_BG_H1_2005 = np.array([0.20, 0.19, 0.20, 0.24, 0.20, 0.21, 0.22, 0.16, 0.19, 0.14, 0.20])
color_H1_2005='red'

# Data H1 2013
W_H1_2013 = np.array([55, 78])
BG_H1_2013 = np.array([4.3, 4.88])
e_BG_H1_2013 = np.array([0.2, 0.15])
color_H1_2013='orange'

# Data ZEUS elec
W_ZEUS_elec = np.array(   [35, 60, 80, 110, 140, 200, 260])
BG_ZEUS_elec = np.array(   [3.55, 3.86,3.98,4.48,4.30,4.65,4.05])
e_BG_ZEUS_elec = np.array([0.27, 0.18, 0.23, 0.22, 0.19, 0.20, 0.38])
color_ZEUS_elec='blue'

# Data ZEUS muon
W_ZEUS_muon = np.array([40,60,80,100,120,140,160])
BG_ZEUS_muon = np.array([3.93,4.02,4.27,4.22,4.28,4.46,4.58])
e_BG_ZEUS_muon = np.array([0.12,0.15,0.15,0.17,0.20,0.27,0.41])
color_ZEUS_muon='green'

fig = plt.figure(figsize=(8, 6))
plt.errorbar(W_H1_2005, BG_H1_2005, yerr=e_BG_H1_2005,  marker='o', linestyle='', color=color_H1_2005, capsize=5, label='H1 (2005)')
plt.errorbar(W_H1_2013, BG_H1_2013, yerr=e_BG_H1_2013,  marker='o', linestyle='', color=color_H1_2013, capsize=5, label='H1 (2013)')
plt.errorbar(W_ZEUS_elec, BG_ZEUS_elec, yerr=e_BG_ZEUS_elec,  marker='o', linestyle='', color=color_ZEUS_elec, capsize=5, label='ZEUS ($e^{-}e^{+}$) (2002)')
plt.errorbar(W_ZEUS_muon, BG_ZEUS_muon, yerr=e_BG_ZEUS_muon,  marker='o', linestyle='', color=color_ZEUS_muon, capsize=5, label=r'ZEUS ($\mu^{-}\mu^{+}$) (2002)')
plt.errorbar(W_values_CLAS12, BG_CLAS12, yerr=BG_errors_CLAS12,  marker='o', linestyle='', color=color_CLAS12, capsize=5, label='CLAS12 (this work)')

#plt.errorbar(mA_values_GlueX, Eg_values_GlueX, xerr=mA_errors_GlueX, yerr=[eg_errors_down_GlueX, eg_errors_up_GlueX], marker='o', linestyle='', color=color_glueX, capsize=5, label='GlueX (2023)')
#plt.errorbar(mA_values_CLAS12, Eg_values_CLAS12, xerr=mA_errors_CLAS12, yerr=[eg_errors_down_CLAS12, eg_errors_up_CLAS12], marker='o', linestyle='', color=color_CLAS12, capsize=5, label='CLAS12')
plt.xlabel('W [GeV]')
plt.ylabel('B$_G$ [GeV$^{-2}$]')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().set_xscale('log')
plt.gca().set_ylim(top=6.)
plt.gca().set_xlim(right=300)
plt.legend(frameon=False, prop={'size': 15})
#plt.legend(loc='upper left', ncol=2,  frameon=False, prop={'size': 15})
plt.text(0.5, 0.5, 'Preliminary', fontsize=80, rotation=45, ha='center', va='center', alpha=0.2, transform=fig.transFigure)


# Save the figure as a PNG file
plt.savefig('Figure_BG.png', dpi=300)


#######################

# Create a square grid with two columns for the two plots
fig = plt.figure(figsize=(12, 8))
gs = GridSpec(1, 2, width_ratios=[1, 1])

# Plot mA as a function of Eg with rotated graph
ax0 = plt.subplot(gs[0])
ax0.errorbar(BG_H1_2005, w_to_tau(W_H1_2005), xerr=e_BG_H1_2005,  marker='o', linestyle='', color=color_H1_2005, capsize=5, label='H1 (2005)')
ax0.errorbar(BG_H1_2013, w_to_tau(W_H1_2013), xerr=e_BG_H1_2013,  marker='o', linestyle='', color=color_H1_2013, capsize=5, label='H1 (2013)')
ax0.errorbar(BG_ZEUS_elec, w_to_tau(W_ZEUS_elec), xerr=e_BG_ZEUS_elec,  marker='o', linestyle='', color=color_ZEUS_elec, capsize=5, label='ZEUS ($e^{-}e^{+}$) (2002)')
ax0.errorbar(BG_ZEUS_muon, w_to_tau(W_ZEUS_muon), xerr=e_BG_ZEUS_muon,  marker='o', linestyle='', color=color_ZEUS_muon, capsize=5, label=r'ZEUS ($\mu^{-}\mu^{+}$) (2002)')
ax0.errorbar(BG_CLAS12, w_to_tau(W_values_CLAS12), xerr=BG_errors_CLAS12,  marker='o', linestyle='', color=color_CLAS12, capsize=5, label='CLAS12 (this work)')
#ax0.errorbar(BG_CLAS12, w_to_tau(W_values_CLAS12), xerr=BG_errors_CLAS12,  marker='', alpha = 0.0, linestyle='', color=color_CLAS12, capsize=5, label='')
ax0.set_xlabel('B$_G$ (GeV$^{-2}$)', fontsize=16)  
ax0.set_ylabel('$x_B$', fontsize=16)  

# Plot mass radius as a function of Eg with rotated graph
ax1 = plt.subplot(gs[1], sharey=ax0)
ax1.errorbar(convert_Gev_to_fm*np.sqrt(2*BG_H1_2005), w_to_tau(W_H1_2005), xerr=e_BG_H1_2005*convert_Gev_to_fm/(2*np.sqrt(2*BG_H1_2005)),  marker='o', linestyle='', color=color_H1_2005, capsize=5, label='H1 (2005)')
ax1.errorbar(convert_Gev_to_fm*np.sqrt(2*BG_H1_2013), w_to_tau(W_H1_2013), xerr=e_BG_H1_2013*convert_Gev_to_fm/(2*np.sqrt(2*BG_H1_2013)),  marker='o', linestyle='', color=color_H1_2013, capsize=5, label='H1 (2013)')
ax1.errorbar(convert_Gev_to_fm*np.sqrt(2*BG_ZEUS_elec), w_to_tau(W_ZEUS_elec), xerr=e_BG_ZEUS_elec*convert_Gev_to_fm/(2*np.sqrt(2*BG_ZEUS_elec)),  marker='o', linestyle='', color=color_ZEUS_elec, capsize=5, label='ZEUS ($e^{-}e^{+}$) (2002)')
ax1.errorbar(convert_Gev_to_fm*np.sqrt(2*BG_ZEUS_muon), w_to_tau(W_ZEUS_muon), xerr=e_BG_ZEUS_muon*convert_Gev_to_fm/(2*np.sqrt(2*BG_ZEUS_muon)),  marker='o', linestyle='', color=color_ZEUS_muon, capsize=5, label=r'ZEUS ($\mu^{-}\mu^{+}$) (2002)')
ax1.errorbar(convert_Gev_to_fm*np.sqrt(2*BG_CLAS12), w_to_tau(W_values_CLAS12), xerr=BG_errors_CLAS12*convert_Gev_to_fm/(2*np.sqrt(2*BG_CLAS12)),  marker='o', linestyle='', color=color_CLAS12, capsize=5, label='CLAS12 (this work)')
#ax1.errorbar(convert_Gev_to_fm*np.sqrt(2*BG_CLAS12), w_to_tau(W_values_CLAS12), xerr=BG_errors_CLAS12*convert_Gev_to_fm/(2*np.sqrt(2*BG_CLAS12)),  alpha = 0.0, marker='', linestyle='', color=color_CLAS12, capsize=5, label='')
ax1.set_xlabel('Gluonic transverse radius (fm)', fontsize=16)
ax1.tick_params(axis='y', labelleft=False)  # Remove y-axis label on the right plot

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
ax0.set_yscale('log')
ax0.set_ylim(ax0.get_ylim()[0], 1)

# Force the x-axis to finish at 0.7
ax1.set_xlim(ax1.get_xlim()[0], 0.9)
ax1.axvline(x=0.84, color='blue', linewidth=2)
ax1.text(0.83, 0.007, 'CODATA charge radius of the proton', color='blue', fontsize=12, rotation=90, va='center', ha='right')





plt.text(0.5, 0.5, 'Preliminary', fontsize=100, rotation=45, ha='center', va='center', alpha=0.2, transform=fig.transFigure)

# Adjust layout to prevent overlap and remove space between graphs
plt.subplots_adjust(wspace=0)

# Save the figure as a PNG file
plt.savefig('BG_and_radius.png', dpi=500)


# Function definition
def r2Pr(r, C0, m):
    return (1/(6*Mn)) * (( 4 * C0 * m**3) / (32*np.pi)) * (r**2) * (m**2) * (m*r - 3) * np.exp(-m*r)

#def min_max_r2Pr(m_min, m_max, C0_min, C0_max):
#    # Constants
#    Mn = 1  # Example value, replace as needed
#    r_values = np.linspace(0., 2, 100)/convert_Gev_to_fm
#    
#
#    ## Function definition
#    #def r2Pr(r, m, C0):
#    #    return (1/(6*Mn)) * (( 4 * C0 * m**3) / (32*np.pi)) * (r**2) * (m**2) * (m*r - 3) * np.exp(-m*r)
#
#    # Initialize lists for storing min and max values
#    min_r2Pr_values = []
#    max_r2Pr_values = []
#
#    # Loop over r values
#    for r in r_values:
#        # Define ranges for m and C0
#        m_vals = np.linspace(m_min, m_max, 100)
#        C0_vals = np.linspace(C0_min, C0_max, 100)
#
#        # Evaluate r2Pr over the grid of (m, C0)
#        values = np.array([r2Pr(r, m, C0) for m in m_vals for C0 in C0_vals])
#        min_r2Pr_values.append(np.min(values))
#        max_r2Pr_values.append(np.max(values))
#
#    return r_values*convert_Gev_to_fm, min_r2Pr_values, max_r2Pr_values

def r2Sr(r, C0, m):
        return (-1/(4*Mn)) * ((4 * C0 * m**3) / (32*np.pi)) * (r**2) * (m**3) * r * np.exp(-m*r)

#def min_max_r2Sr(m_min, m_max, C0_min, C0_max):
#    # Constants
#    Mn = 1  # Example value, replace as needed
#    r_values = np.linspace(0., 2, 100)/convert_Gev_to_fm
#    
#
#    ## Function definition
#    #def r2Sr(r, m, C0):
#    #    return (-1/(4*Mn)) * ((4 * C0 * m**3) / (32*np.pi)) * (r**2) * (m**3) * r * np.exp(-m*r)
#
#    # Initialize lists for storing min and max values
#    min_r2Sr_values = []
#    max_r2Sr_values = []
#
#    # Loop over r values
#    for r in r_values:
#        # Define ranges for m and C0
#        m_vals = np.linspace(m_min, m_max, 100)
#        C0_vals = np.linspace(C0_min, C0_max, 100)
#
#        # Evaluate r2Pr over the grid of (m, C0)
#        values = np.array([r2Sr(r, m, C0) for m in m_vals for C0 in C0_vals])
#        min_r2Sr_values.append(np.min(values))
#        max_r2Sr_values.append(np.max(values))
#
#    return r_values*convert_Gev_to_fm, min_r2Sr_values, max_r2Sr_values

##### Presure and shear force disribution #####

def Plot_PandS(GPD=False):
    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(1, 2, width_ratios=[1, 1])
    if(GPD):
        min_r = 0.0
        max_r = 4.0

    else:
        min_r = 0.0
        max_r = 2.0

    r_values=np.linspace(min_r, max_r, 100)
    r_GeV = r_values/convert_Gev_to_fm

    if(fix_mC):
        corr_GPD_B = 0.0
        corr_Holo_B = 0.0

    if(GPD):
        #corr = 0.0
        #if(not fix_C0  and not fix_mA):
        #    corr=result_GPD_C_B.covar[0][2]/np.sqrt(result_GPD_C_B.covar[0][0]*result_GPD_C_B.covar[2][2])
        #elif (fix_mA):
        #    corr=result_GPD_C_B.covar[0][1]/np.sqrt(result_GPD_C_B.covar[1][1]*result_GPD_C_B.covar[0][0])
        min_r2Pr_values_B_C, max_r2Pr_values_B_C = compute_uncertainty_bands(r2Pr, x0 = result_GPD_C_B.params['C_0'].value, y0 = 0.91, sigma_x=result_GPD_C_B.params['C_0'].stderr, sigma_y=0.10, rho=corr_GPD_B, z_vals=r_GeV)

        min_r2Sr_values_B_C, max_r2Sr_values_B_C = compute_uncertainty_bands(r2Sr, x0 = result_GPD_C_B.params['C_0'].value, y0 = 0.91, sigma_x=result_GPD_C_B.params['C_0'].stderr, sigma_y=0.10, rho=corr_GPD_B, z_vals=r_GeV)
        
        min_r2Pr_values_all, max_r2Pr_values_all = compute_uncertainty_bands(r2Pr, x0 = result_GPD.params['C_0'].value, y0 =result_GPD.params['m_C'].value, sigma_x=result_GPD.params['C_0'].stderr, sigma_y=result_GPD.params['m_C'].stderr, rho=result_GPD.covar[0][2]/np.sqrt(result_GPD.covar[0][0]*result_GPD.covar[2][2]), z_vals=r_GeV)
        
        min_r2Sr_values_all, max_r2Sr_values_all = compute_uncertainty_bands(r2Sr, x0 = result_GPD.params['C_0'].value, y0 =result_GPD.params['m_C'].value, sigma_x=result_GPD.params['C_0'].stderr, sigma_y=result_GPD.params['m_C'].stderr, rho=result_GPD.covar[0][2]/np.sqrt(result_GPD.covar[0][0]*result_GPD.covar[2][2]), z_vals=r_GeV)
        
        min_r2Pr_values_C, max_r2Pr_values_C = compute_uncertainty_bands(r2Pr, x0 = -1.478, y0 =0.774, sigma_x=0.467, sigma_y=0.092, rho=1., z_vals=r_GeV)

        min_r2Sr_values_C, max_r2Sr_values_C = compute_uncertainty_bands(r2Sr,  x0 = -1.478, y0 =0.774, sigma_x=0.467, sigma_y=0.092, rho=1., z_vals=r_GeV)

        

        ax0 = plt.subplot(gs[0])
        ax0.plot([], [], ' ', label=r'$\mathbf{GPD\ Model}$') 
        #ax0.fill_between(t_values, A_g_lower_HC_GPD, A_g_upper_HC_GPD, color=color_hall_c, alpha=0.5, label='Hall C only (Nature)')
        ax0.fill_between(r_values, min_r2Pr_values_B_C/convert_Gev_to_fm, max_r2Pr_values_B_C/convert_Gev_to_fm, color='purple', alpha=0.5, label=r'Hall B ($A(0)$ and $m_C$ fixed)')
        ax0.fill_between(r_values, min_r2Pr_values_C/convert_Gev_to_fm, max_r2Pr_values_C/convert_Gev_to_fm, color=color_hall_c, alpha=0.5, label=r'Hall C ($A(0)$ fixed)')
        ax0.fill_between(r_values, min_r2Pr_values_all/convert_Gev_to_fm, max_r2Pr_values_all/convert_Gev_to_fm, color='blue', alpha=0.5, label=r'All data combined ($A(0)$ fixed)')

        ax0.axhline(0, color="black", linestyle="--", linewidth=1) 


        ax1 = plt.subplot(gs[1])
        #ax1.fill_between(t_values, D_g_lower_HC_GPD, D_g_upper_HC_GPD, color=color_hall_c, alpha=0.5, label='Hall C only')
        ax1.fill_between(r_values, min_r2Sr_values_B_C/convert_Gev_to_fm, max_r2Sr_values_B_C/convert_Gev_to_fm, color='purple', alpha=0.5, label=r'Hall B ($A(0)$ and $m_C$ fixed)')
        ax1.fill_between(r_values, min_r2Sr_values_C/convert_Gev_to_fm, max_r2Sr_values_C/convert_Gev_to_fm, color=color_hall_c, alpha=0.5, label=r'Hall C ($A(0)$ fixed)')
        ax1.fill_between(r_values, min_r2Sr_values_all/convert_Gev_to_fm, max_r2Sr_values_all/convert_Gev_to_fm, color='blue', alpha=0.5, label=r'All data combined ($A(0)$ fixed)')

        ax1.axhline(0, color="black", linestyle="--", linewidth=1) 

    else:       
                                                                   
        C_0_C = -0.45 
        C_0_error_C = 0.132
        m_C_C = 1.12 
        m_C_error_C = 0.21
    
        min_r2Pr_values_B_C, max_r2Pr_values_B_C = compute_uncertainty_bands(r2Pr, x0 = result_Holo_C_B.params['C_0'].value, y0 = 1.12, sigma_x=result_Holo_C_B.params['C_0'].stderr, sigma_y=0.21, rho=corr_Holo_B, z_vals=r_GeV)

        min_r2Sr_values_B_C, max_r2Sr_values_B_C = compute_uncertainty_bands(r2Sr, x0 = result_Holo_C_B.params['C_0'].value, y0 = 1.12, sigma_x=result_Holo_C_B.params['C_0'].stderr, sigma_y=0.21, rho=corr_Holo_B, z_vals=r_GeV)
        
        min_r2Pr_values_all, max_r2Pr_values_all = compute_uncertainty_bands(r2Pr, x0 = result_holo.params['C_0'].value, y0 =result_holo.params['m_C'].value, sigma_x=result_holo.params['C_0'].stderr, sigma_y=result_holo.params['m_C'].stderr, rho=result_holo.covar[0][2]/np.sqrt(result_holo.covar[0][0]*result_holo.covar[2][2]), z_vals=r_GeV)
        
        min_r2Sr_values_all, max_r2Sr_values_all = compute_uncertainty_bands(r2Sr, x0 = result_holo.params['C_0'].value, y0 =result_holo.params['m_C'].value, sigma_x=result_holo.params['C_0'].stderr, sigma_y=result_holo.params['m_C'].stderr, rho=result_holo.covar[0][2]/np.sqrt(result_holo.covar[0][0]*result_holo.covar[2][2]), z_vals=r_GeV)
        
        min_r2Pr_values_C, max_r2Pr_values_C = compute_uncertainty_bands(r2Pr, x0 = C_0_C, y0 =m_C_C, sigma_x=C_0_error_C, sigma_y=m_C_error_C, rho=0.9, z_vals=r_GeV)

        #min_r2Pr_values_C, max_r2Pr_values_C = compute_uncertainty_bands(r2Pr, x0 = -1.93/4., y0 =1.07, sigma_x=0, sigma_y=0, rho=0.0, z_vals=r_GeV)

        min_r2Sr_values_C, max_r2Sr_values_C = compute_uncertainty_bands(r2Sr, x0 = C_0_C, y0 =m_C_C, sigma_x=C_0_error_C, sigma_y=m_C_error_C, rho=0.9, z_vals=r_GeV)

        

        ax0 = plt.subplot(gs[0])
        ax0.plot([], [], ' ', label=r'$\mathbf{Holographic\ Model}$') 
        #ax0.fill_between(t_values, A_g_lower_HC_GPD, A_g_upper_HC_GPD, color=color_hall_c, alpha=0.5, label='Hall C only (Nature)')
        ax0.fill_between(r_values, min_r2Pr_values_B_C/convert_Gev_to_fm, max_r2Pr_values_B_C/convert_Gev_to_fm, color='purple', alpha=0.5, label=r'Hall B ($A(0)$ and $m_C$ fixed)')
        ax0.fill_between(r_values, min_r2Pr_values_C/convert_Gev_to_fm, max_r2Pr_values_C/convert_Gev_to_fm, color=color_hall_c, alpha=0.5, label=r'Hall C ($A(0)$ fixed)')
        ax0.fill_between(r_values, min_r2Pr_values_all/convert_Gev_to_fm, max_r2Pr_values_all/convert_Gev_to_fm, color='blue', alpha=0.5, label=r'All data combined ($A(0)$ fixed)')

        ax0.axhline(0, color="black", linestyle="--", linewidth=1) 


        ax1 = plt.subplot(gs[1])
        #ax1.fill_between(t_values, D_g_lower_HC_GPD, D_g_upper_HC_GPD, color=color_hall_c, alpha=0.5, label='Hall C only')
        ax1.fill_between(r_values, min_r2Sr_values_B_C/convert_Gev_to_fm, max_r2Sr_values_B_C/convert_Gev_to_fm, color='purple', alpha=0.5, label=r'Hall B ($A(0)$ and $m_C$ fixed)')
        ax1.fill_between(r_values, min_r2Sr_values_C/convert_Gev_to_fm, max_r2Sr_values_C/convert_Gev_to_fm, color=color_hall_c, alpha=0.5, label=r'Hall C ($A(0)$ fixed)')
        ax1.fill_between(r_values, min_r2Sr_values_all/convert_Gev_to_fm, max_r2Sr_values_all/convert_Gev_to_fm, color='blue', alpha=0.5, label=r'All data combined ($A(0)$ fixed)')

        ax1.axhline(0, color="black", linestyle="--", linewidth=1) 


   
    ax0.set_xlabel('r [fm]')
    ax0.set_ylabel(r'$r^{2}p_{g}(r)$ [GeV/fm]')
    ax0.set_xlim(left=min_r)
    ax0.set_xlim(right=max_r)
    #ax0.set_yscale('log')


    ax1.set_xlabel('r [fm]')
    ax1.set_ylabel(r'$r^{2}s_{g}(r)$ [GeV/fm]')
    ax1.set_xlim(left=min_r)
    ax1.set_xlim(right=max_r)



    ax0.legend(frameon=False)
    #ax1.legend(frameon=False)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    if(GPD):
        fig.savefig("PandS_GPD.png", dpi=500)
    else:
        fig.savefig("PandS_Holo.png", dpi=500)
    #plt.show()

Plot_PandS(GPD=False)
Plot_PandS(GPD=True)