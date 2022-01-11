import math
import pandas as pd
import sys

sys.path.insert(0, '/mnt/c/Users/Lukas/cernbox/Data/Analysis_Code/bin/')

Masses = pd.read_fwf('/mnt/c/Users/Lukas/cernbox/Data/Analysis_Code/bin/mass20.txt', #usecols=(2,3,4,6,9,10,11,12,18,21,20),
              names=['1', 'N-Z', 'N', 'Z', 'A', 'Unnamed: 5', 'element', 'O', 'Unnamed: 8',
                       'mass_excess', 'mass_excess_err', 'ebinding', 'nan1', 'ebinding_err', 'nan2', 'ET',
                       'beta_decay_energy', 'beta_decay_energy_err', 'nan18', 'atomic_mass_raw', 'nan20',
                       'atomic_mass_comma', 'atomic_mass_err'],
              widths=(1,3,5,5,5,1,3,4,1,14,12,13,1,10,1,2,13,11,1,3,1,13,12),
              header=28,
              index_col=False)
#
c = 299792458 # m/s
e = 1.602176634e-19 
u = 931494.10242 # keV/c**2
u_err = 000.00000028 # MeV/c**2

m_p = float(Masses["atomic_mass_raw"][(Masses["A"]==1) & (Masses["element"]=="H")]) + float(Masses["atomic_mass_comma"][(Masses["A"]==1) & (Masses["element"]=="H")])/1e6
m_n = float(Masses["atomic_mass_raw"][(Masses["A"]==1) & (Masses["element"]=="n")]) + float(Masses["atomic_mass_comma"][(Masses["A"]==1) & (Masses["element"]=="n")])/1e6
m_85_Rb = float(Masses["atomic_mass_raw"][(Masses["A"]==85) & (Masses["element"]=="Rb")]) + float(Masses["atomic_mass_comma"][(Masses["A"]==85) & (Masses["element"]=="Rb")])/1e6
m_85_Rb_err = float(Masses["atomic_mass_err"][(Masses["A"]==85) & (Masses["element"]=="Rb")])
m_133_Cs = float(Masses["atomic_mass_raw"][(Masses["A"]==133) & (Masses["element"]=="Cs")]) + float(Masses["atomic_mass_comma"][(Masses["A"]==133) & (Masses["element"]=="Cs")])/1e6
m_133_Cs_err = float(Masses["atomic_mass_err"][(Masses["A"]==133) & (Masses["element"]=="Cs")])


def calc_weighted_average(x, s):
    """ Takes array of values (x) plus array of errors (s) and returns weighted average"""
    if len(x) != len(s):
        print("Arrays must have same length")
        return 0
    # 
    sum_avg = 0
    sum_w = 0
    for i in range(len(x)):
        w = s[i]**(-2)
        sum_w += w
        sum_avg += w * x[i]
    #
    return sum_avg/sum_w
def calc_weighted_averge_err(s):
    """Takes array s of individual errors"""
    sum_w = 0
    for i in range(len(s)):
        sum_w += s[i]**(-2)
    #
    return math.sqrt(1/sum_w)

def calc_red_chi_square(x, s):
    """Calculates reduced chi square for array of values (x) and array of errors (s)"""
    if len(x) != len(s):
        print("Arrays must have same length")
        return 0
    # 
    weighted_average = calc_weighted_average(x, s)
    #
    chi_square = 0
    for i in range(len(x)):
        chi_square += (x[i]-weighted_average)**2 / s[i]**2
    #
    return chi_square / (len(x)-1)
    

def calc_C_ToF(t,t1,t2):
    """ Calculation of C_ToF based on Franks 2013 Nature paper"""
    return (2*t-t1-t2)/(2*(t1-t2))

def calc_C_ToF_err(t, t_err, t1, t1_err, t2, t2_err):
    """"""
    del_t = 1 / (t1-t2)
    del_t1 = (-t+t2)/(t1-t2)**2
    del_t2 = (t-t1)/(t1-t2)**2
    #
    return math.sqrt( (del_t*t_err)**2 + (del_t1*t1_err)**2 + (del_t2*t2_err)**2 )
    
def calc_sqrt_m(C_tof, m1, m2):
    """ Calculation of the sqrt of the mass of the species of interest"""
    Sigma_ref = math.sqrt(m1) + math.sqrt(m2)
    Delta_ref = math.sqrt(m1) - math.sqrt(m2)
    #
    return C_tof * Delta_ref + Sigma_ref/2 

def calc_sqrt_m_err(C_tof, C_tof_err, m1, m1_err, m2, m2_err):
    """ Calculation of the err on the sqrt of the mass"""
    del_C_tof = math.sqrt(m1) - math.sqrt(m2)
    del_m1 = C_tof + 1/2
    del_m2 = - C_tof + 1/2
    #
    sqrt_m1_err = 1/2 * m1**(-1/2) * m1_err
    sqrt_m2_err = 1/2 * m2**(-1/2) * m2_err
    #
    return math.sqrt( (del_C_tof * C_tof_err)**2 + ( del_m1 * sqrt_m1_err )**2 + ( del_m2 * sqrt_m2_err )**2 )

def calc_m_err_alternative(sqrt_m, sqrt_m_err):
    """ Calculation of the mass error using the error on the sqrt of the mass"""
    return 2 * sqrt_m * sqrt_m_err

def calc_m_err(C_tof, C_tof_err, m1, m1_err, m2, m2_err):
    """ Calculation of the mass error using error propagation"""
    delta_ref = math.sqrt(m1) - math.sqrt(m2)
    sigma_ref = math.sqrt(m1) + math.sqrt(m2)
    #
    del_C_tof = 2 * C_tof * delta_ref**2 + delta_ref * sigma_ref
    del_m1 = C_tof**2 * (1-m1**(-1/2)) + C_tof + 1/4 * (1+m1**(-1/2))
    del_m2 = C_tof**2 * (1-m2**(-1/2)) + C_tof + 1/4 * (1+m2**(-1/2))
    #
    return math.sqrt( (del_C_tof*C_tof_err)**2 + (del_m1 * m1_err)**2 + (del_m2 * m2_err)**2 )