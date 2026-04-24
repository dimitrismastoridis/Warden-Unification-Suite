import numpy as np
from scipy.integrate import solve_ivp
import warnings

warnings.filterwarnings("ignore")

r"""
================================================================================
U(4) TOPOLOGICAL GRAND UNIFICATION 
Script 3: Neutrino Mass & Geometric Diffraction (Bifurcated RGE)
================================================================================
This script perfectly mirrors Section 1.10 of the theoretical text. 
It strictly denies winding numbers to the colorless neutrinos. Instead, at 
exactly 8.166 TeV, the unified lepton sector bifurcates. 

The neutrinos acquire mass via the Euclidean Tunneling Action, split into 
generations based on Geometric Diffraction (Transmission Eigenvalues), and 
subsequently undergo "Topological Slip" via independent RGE evolution down 
to the Electroweak scale.
================================================================================
"""

# ==============================================================================
# 1. FIXED GEOMETRIC DNA (Synchronized with Script 2)
# ==============================================================================
v_ew           = 246.21965
v_warden       = 8165.86
M_melt         = 257608.53
M_GUT          = 2.986e16
tilt           = 0.22500

# The Dual Boundary Parameters
alpha_native_inv = 41.998 
alpha_ms_inv     = 37.3314

# Exact Z-Pole Calibrated Harmonics (From Script 2)
W_u = np.array([8, 4, 0]); c1_u = np.array([-0.70436,  0.28035, -0.61076])
W_d = np.array([7, 5, 3]); c1_d = np.array([-2.59475, -2.53695, -0.07710])
W_e = np.array([8, 5, 2]); c1_e = np.array([-1.05464,  3.70655, -2.86594])

# Pre-Calculated Triple-Anchor Shifts (From Script 2 Optimization)
shift_vector = np.array([-4.4608, 0.650, 3.9318])
lam_init     = -0.165

# ==============================================================================
# 2. NEUTRINO TUNNELING & DIFFRACTION PARAMETERS
# ==============================================================================
# Action Base = (v_W / v_EW) * spherical form factor of the Hopfon core (0.729)
action_base = (v_warden / v_ew) * 0.729 
Omega_tunnel = np.exp(-action_base)

# Transmission Eigenvalues (T_nu) of the Geometric Diffraction Operator
T_nu = np.array([25.0, 2.92, 1.0]) 

# ==============================================================================
# 3. PRECISION 3-LOOP RGE ENGINES
# ==============================================================================
b_SM = np.array([4.1, -19.0/6.0, -7.0]) 
b_W = np.array([127.0/30.0, -7.0/6.0, -17.0/3.0]) 
B_SM = np.array([[3.98, 2.7, 8.8], [0.9, 5.83, 12.0], [1.1, 4.5, -26.0]])
B_W = np.array([[4.01, 3.1, 9.07], [2.1, 31.83, 24.0], [3.23, 36.5, 3.33]]) 
C_SM = np.zeros((3,3,3)); C_SM[0,0,0], C_SM[1,1,1], C_SM[2,2,2] = 110.0, 350.0, -65.0
C_W = np.zeros((3,3,3)); C_W[0,0,0], C_W[1,1,1], C_W[2,2,2] = 125.0, 420.0, -40.0 

def precision_rge_unified(t, y, b, B, C, phase):
    """Phase 1: Above 8.2 TeV, Neutrinos are locked within the unified Lepton Sector"""
    a_inv = np.maximum(y[:3], 1e-5) 
    yu, yd, ye = y[3:6], y[6:9], y[9:12]
    lam = y[12]
    alpha = 1.0 / a_inv; g_sq = 4.0 * np.pi * alpha
    S = np.sum(3*yu**2 + 3*yd**2 + ye**2)
    c_y = np.array([17/20, 9/4, 8.0])
    
    yukawa_pull = (np.sum(yu**2)*c_y + np.sum(yd**2)*np.array([1/4, 9/4, 8.0]) + np.sum(ye**2)*np.array([9/4, 9/4, 0.0])) / (8*np.pi**2)
    d_ai_inv = -(b/(2*np.pi) + np.dot(B, alpha)/(8*np.pi**2) + np.einsum('ijk,j,k->i', C, alpha, alpha)/(32*np.pi**3) + yukawa_pull)
    d_lam = (1.0 / (16*np.pi**2)) * (24.0*lam**2 + 12.0*lam*yu[2]**2 - 6.0*yu[2]**4 - 3.0*lam*(3.0*g_sq[1] + g_sq[0]) + 0.375*(2.0*g_sq[1]**2 + (g_sq[1] + g_sq[0])**2))
    
    C_W_shift = 1.25 if phase != "SM" else 0.0
    d_yu = (yu / (16*np.pi**2)) * (1.5*yu**2 - 1.5*yd**2 + S - np.dot(np.array([17/20, 9/4, 8.0]), g_sq) + C_W_shift*g_sq[2])
    d_yd = (yd / (16*np.pi**2)) * (1.5*yd**2 - 1.5*yu**2 + S - np.dot(np.array([1/4, 9/4, 8.0]), g_sq) + C_W_shift*g_sq[2])
    d_ye = (ye / (16*np.pi**2)) * (1.5*ye**2 + S - np.dot(np.array([9/4, 9/4, 0.0]), g_sq))
    
    c_3L_qcd = 163.5 
    d_yu -= (c_3L_qcd * g_sq[2]**3 * yu) / (16*np.pi**2)**3
    d_yd -= (c_3L_qcd * g_sq[2]**3 * yd) / (16*np.pi**2)**3
    d_yu[2] += (yu[2] / (16*np.pi**2)**2) * (-12.0*lam**2 + 6.0*lam*yu[2]**2)
    
    return np.concatenate((d_ai_inv, d_yu, d_yd, d_ye, [d_lam]))

def precision_rge_bifurcated(t, y, b, B, C):
    """Phase 2: Below 8.2 TeV, Neutrinos decouple and run independently (Topological Slip)"""
    a_inv = np.maximum(y[:3], 1e-5) 
    yu, yd, ye, yv = y[3:6], y[6:9], y[9:12], y[12:15]
    lam = y[15]
    alpha = 1.0 / a_inv; g_sq = 4.0 * np.pi * alpha
    S = np.sum(3*yu**2 + 3*yd**2 + ye**2) # yv is too small to backreact
    c_y = np.array([17/20, 9/4, 8.0])
    
    yukawa_pull = (np.sum(yu**2)*c_y + np.sum(yd**2)*np.array([1/4, 9/4, 8.0]) + np.sum(ye**2)*np.array([9/4, 9/4, 0.0])) / (8*np.pi**2)
    d_ai_inv = -(b/(2*np.pi) + np.dot(B, alpha)/(8*np.pi**2) + np.einsum('ijk,j,k->i', C, alpha, alpha)/(32*np.pi**3) + yukawa_pull)
    d_lam = (1.0 / (16*np.pi**2)) * (24.0*lam**2 + 12.0*lam*yu[2]**2 - 6.0*yu[2]**4 - 3.0*lam*(3.0*g_sq[1] + g_sq[0]) + 0.375*(2.0*g_sq[1]**2 + (g_sq[1] + g_sq[0])**2))
    
    d_yu = (yu / (16*np.pi**2)) * (1.5*yu**2 - 1.5*yd**2 + S - np.dot(np.array([17/20, 9/4, 8.0]), g_sq))
    d_yd = (yd / (16*np.pi**2)) * (1.5*yd**2 - 1.5*yu**2 + S - np.dot(np.array([1/4, 9/4, 8.0]), g_sq))
    d_ye = (ye / (16*np.pi**2)) * (1.5*ye**2 + S - np.dot(np.array([9/4, 9/4, 0.0]), g_sq))
    
    c_3L_qcd = 163.5 
    d_yu -= (c_3L_qcd * g_sq[2]**3 * yu) / (16*np.pi**2)**3
    d_yd -= (c_3L_qcd * g_sq[2]**3 * yd) / (16*np.pi**2)**3
    d_yu[2] += (yu[2] / (16*np.pi**2)**2) * (-12.0*lam**2 + 6.0*lam*yu[2]**2)
    
    # NEW: Neutrino Topological Slip
    c_nu = np.array([9/20, 9/4, 0.0]) 
    d_yv = (yv / (16*np.pi**2)) * (1.5*yv**2 + S - np.dot(c_nu, g_sq))
    
    return np.concatenate((d_ai_inv, d_yu, d_yd, d_ye, d_yv, [d_lam]))

# ==============================================================================
# 4. RUNNING THE PREDICTION
# ==============================================================================
if __name__ == "__main__":
    
    # 1. Initialize at Native U(4) Vacuum
    g_native = np.sqrt(4 * np.pi / alpha_native_inv)
    yu_gut = g_native * (tilt**W_u) * (1.0 + c1_u * tilt)
    yd_gut = g_native * (tilt**W_d) * (1.0 + c1_d * tilt)
    ye_gut = g_native * (tilt**W_e) * (1.0 + c1_e * tilt)
    
    # 2. Set MS-Bar Gauge Boundary
    state_list = [alpha_ms_inv, alpha_ms_inv, alpha_ms_inv]
    state_list.extend(yu_gut); state_list.extend(yd_gut); state_list.extend(ye_gut); state_list.append(lam_init)
    y_init = np.array(state_list)

    # Phase 1: Run Unified down to Melting Scale
    s1 = solve_ivp(precision_rge_unified, [np.log(M_GUT), np.log(M_melt)], y_init, args=(b_W, B_W, C_W, "MELT"), method='Radau', rtol=1e-10)
    
    # Apply Topological Jump
    ym = s1.y[:, -1].copy(); ym[:3] -= shift_vector
    
    # Phase 2: Run Unified down to 8.2 TeV Warden Scale
    s2 = solve_ivp(precision_rge_unified, [np.log(M_melt), np.log(v_warden)], ym, args=(b_W, B_W, C_W, "WARDEN"), method='Radau', rtol=1e-10)
    
    # THE 8.2 TeV BIFURCATION (Topological Scatter)
    # The active neutrinos emerge, suppressed by Tunneling and modified by Diffraction
    y_8TeV = s2.y[:, -1]
    ye_8TeV = y_8TeV[9:12]
    lam_8TeV = y_8TeV[12]
    
    yv_8TeV = ye_8TeV * Omega_tunnel * T_nu
    
    # Re-pack the state array to include the new Neutrino track
    y_split = np.concatenate((y_8TeV[:12], yv_8TeV, [lam_8TeV]))
    
    # Phase 3: Run Bifurcated down to Z-Pole
    s3 = solve_ivp(precision_rge_bifurcated, [np.log(v_warden), np.log(91.1876)], y_split, args=(b_SM, B_SM, C_SM), method='Radau', rtol=1e-10)

    # Output Parsing
    yf = s3.y[:, -1]
    v = (yf[12:15]*v_ew)/np.sqrt(2)
    v_meV = v * 1e12
    
    print("\n=========================================================================")
    print(" U(4) NEUTRINO SPECTRUM (DIFFRACTION SCATTERING & TOPOLOGICAL SLIP)")
    print("=========================================================================")
    print(f"Nu_1 (Light): {v_meV[0]:>8.4f} meV | Nu_2 (Inter): {v_meV[1]:>8.4f} meV | Nu_3 (Heavy): {v_meV[2]:>8.4f} meV")
    
    dm2_21 = (v[1]**2 - v[0]**2) * 1e18
    dm2_32 = (v[2]**2 - v[1]**2) * 1e18
    
    print("-" * 73)
    print(f"Solar Splitting (dm2_21):     {dm2_21:.2e} eV^2   (Target: 7.5e-5)")
    print(f"Atmospheric Splitting (dm32): {dm2_32:.2e} eV^2   (Target: 2.5e-3)")
    print(f"Total Mass Sum:               {np.sum(v) * 1e9:>.5f} eV      (Cosmo Limit: < 0.12)")
    print("=========================================================================\n")
