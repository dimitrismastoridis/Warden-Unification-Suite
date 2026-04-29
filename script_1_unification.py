import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

warnings.filterwarnings("ignore")

# ==============================================================================
# 1. EXPERIMENTAL NOISE INPUTS
# ==============================================================================
mz_input = (91.1876, 0.0021)
a_em_inv_input = (127.955, 0.015)
v_ew_input = (246.21965, 0.00006)
m_top_input = (172.57, 0.40)
as_mz_input = (0.1179, 0.0009) 

m_warden_input = (8165.86, 23.14) 
m_melt_input = (257608.53, 730.0)

b_SM = np.array([4.1, -19.0/6.0, -7.0])
b_W = np.array([127.0/30.0, -7.0/6.0, -17.0/3.0])
B_SM = np.array([[199.0/50.0, 2.7, 8.8], [0.9, 5.83, 12.0], [1.1, 4.5, -26.0]])
B_W = np.array([[4.01, 3.1, 9.07], [2.1, 31.83, 24.0], [3.23, 36.5, 3.33]])
C_SM = np.zeros((3,3,3)); C_SM[0,0,0], C_SM[1,1,1], C_SM[2,2,2] = 110.0, 350.0, -65.0
C_W = np.zeros((3,3,3)); C_W[0,0,0], C_W[1,1,1], C_W[2,2,2] = 125.0, 420.0, -40.0
c_t = np.array([1.7, 1.5, 2.0])
delta_sm_mz = np.array([0.0012, 0.0045, -0.0110]) 

warden_weights = np.array([2.0/15.0, 2.0, 4.0/3.0]) 

# FROZEN ARRAY: The extracted Cho-Savvidy Imprint
Theta_D = np.array([0.88794, 1.79871, -1.52049])

# ==============================================================================
# 2. RAW RGE ENGINE 
# ==============================================================================
def rge_precision(t, y, b, B, C):
    a_inv = np.maximum(y[:3], 1e-5)
    yt = np.clip(y[3], 0.0, 10.0)
    alpha = 1.0 / a_inv
    g_sq = 4.0 * np.pi * alpha
    
    l1 = b / (2 * np.pi)
    l2 = np.dot(B, alpha) / (8 * np.pi**2)
    l3 = np.einsum('ijk,j,k->i', C, alpha, alpha) / (32 * np.pi**3)
    yukawa_pull = (yt**2 / (8 * np.pi**2)) * c_t
    
    d_yt = (yt / (16 * np.pi**2)) * (4.5 * yt**2 - 8.0*g_sq[2] - 2.25*g_sq[1] - 0.85*g_sq[0])
    return np.concatenate((-(l1 + l2 + l3 - yukawa_pull), [d_yt]))

def evaluate_universe(s2w_guess, mz, a_em_inv_mz, v_ew, m_top_pole, as_mz, m_warden, m_melt):
    as_pi = as_mz / np.pi
    mt_msbar = m_top_pole * (1.0 - (4.0/3.0)*as_pi - 1.0414*(as_pi**2) - 3.3714*(as_pi**3))
    yt_start = np.sqrt(2) * mt_msbar / v_ew

    phase_width = np.log(m_melt / m_warden)
    Theta_Theory = -(warden_weights / (2 * np.pi)) * phase_width + Theta_D 

    y_mz = np.array([
        (3/5)*a_em_inv_mz*(1-s2w_guess) + delta_sm_mz[0], 
        a_em_inv_mz*s2w_guess + delta_sm_mz[1], 
        1/as_mz + delta_sm_mz[2], 
        yt_start
    ])
    
    s1 = solve_ivp(rge_precision, [np.log(mz), np.log(m_warden)], y_mz, args=(b_SM, B_SM, C_SM), method='Radau')
    s2 = solve_ivp(rge_precision, [np.log(m_warden), np.log(m_melt)], s1.y[:, -1], args=(b_W, B_W, C_W), method='Radau')
    
    ym = s2.y[:, -1].copy()
    ym[:3] += Theta_Theory 
    
    def cross(ln_scale):
        s3 = solve_ivp(rge_precision, [np.log(m_melt), ln_scale], ym, args=(b_W, B_W, C_W), method='Radau')
        return s3.y[0, -1] - s3.y[1, -1]
    
    ln_mgut = brentq(cross, np.log(1e15), np.log(1e18))
    s_f = solve_ivp(rge_precision, [np.log(m_melt), ln_mgut], ym, args=(b_W, B_W, C_W), method='Radau')
    return s_f.y[1, -1] - s_f.y[2, -1], ln_mgut

# ==============================================================================
# 3. BLIND EXECUTION
# ==============================================================================
iterations = 500 
s2w_results, mgut_results, ms_bar_results = [], [], []
failed_runs = 0

print(f"Initiating RAW BLIND U(4) Monte Carlo ({iterations} Universes)...")
np.random.seed(42) 

for i in range(iterations):
    _mz = np.random.normal(*mz_input)
    _a_em = np.random.normal(*a_em_inv_input)
    _v_ew = np.random.normal(*v_ew_input)
    _m_top = np.random.normal(*m_top_input)
    _as_mz = np.random.normal(*as_mz_input)
    _m_w = np.random.normal(*m_warden_input)
    _m_m = np.random.normal(*m_melt_input)

    try:
        def objective(x): return evaluate_universe(x, _mz, _a_em, _v_ew, _m_top, _as_mz, _m_w, _m_m)[0]
        opt_s2w = brentq(objective, 0.228, 0.235)
        _, opt_ln_mgut = evaluate_universe(opt_s2w, _mz, _a_em, _v_ew, _m_top, _as_mz, _m_w, _m_m)
        
        s2w_results.append(opt_s2w)
        mgut_results.append(np.exp(opt_ln_mgut))
    except ValueError:
        failed_runs += 1 

print("\n================================================================")
print("     U(4) RAW BLIND PREDICTIONS")
print("================================================================")
print(f"Convergence Rate     : {((iterations - failed_runs)/iterations)*100:.1f}%")
print(f"OUTPUT sin^2 theta_W : {np.mean(s2w_results):.5f} ± {np.std(s2w_results):.5f}")
print(f"OUTPUT M_GUT Scale   : {np.mean(mgut_results):.3e} GeV")
print("================================================================\n")
