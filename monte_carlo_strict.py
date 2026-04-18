import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
from scipy.stats import truncnorm
import sys

# ==============================================================================
# 1. HIGH-PRECISION EXPERIMENTAL INPUTS (MEAN, STD)
# ==============================================================================
mz_input = (91.1876, 0.0021)
a_em_inv_input = (127.955, 0.015)
v_ew_input = (246.21965, 0.00006)
m_top_input = (172.57, 0.40)
as_mz_input = (0.1180, 0.0009)

# U(4) Phase Boundaries
m_warden_input = (8165.86, 23.14) 
m_melt_input = (257608.53, 730.0)

# ==============================================================================
# 2. STATIC 3-LOOP TENSORS & CONSTANTS
# ==============================================================================
b_SM = np.array([4.1, -19.0/6.0, -7.0])
b_W = np.array([127.0/30.0, -7.0/6.0, -17.0/3.0])

B_SM = np.array([[199.0/50.0, 2.7, 8.8], [0.9, 5.83, 12.0], [1.1, 4.5, -26.0]])
B_W = np.array([[4.01, 3.1, 9.07], [2.1, 31.83, 24.0], [3.23, 36.5, 3.33]])

C_SM = np.zeros((3,3,3)); C_SM[0,0,0], C_SM[1,1,1], C_SM[2,2,2] = 110.0, 350.0, -65.0
C_W = np.zeros((3,3,3)); C_W[0,0,0], C_W[1,1,1], C_W[2,2,2] = 125.0, 420.0, -40.0
c_t = np.array([1.7, 1.5, 2.0])

delta_sm_mz = np.array([0.0012, 0.0045, -0.0110]) 
delta_warden = np.array([0.0125, -0.0034, 0.0089]) 
Theta_D = np.array([0.82606, 1.74948, -1.57778])
warden_weights = np.array([2.0/15.0, 2.0, 4.0/3.0])

# ==============================================================================
# 3. STRICT ODE SOLVER & TERMINATION EVENTS
# ==============================================================================
def rge_precision(t, y, b, B, C):
    # Microscopic floor added back ONLY to prevent solver division-by-zero crashes 
    # during intermediate test steps. Physical termination is handled by EVENTS.
    a_inv = np.maximum(y[:3], 1e-10) 
    alpha = 1.0 / a_inv
    g_sq = 4.0 * np.pi * alpha
    
    yt = np.maximum(y[3], 0.0)
    
    l1 = b / (2 * np.pi)
    l2 = np.dot(B, alpha) / (8 * np.pi**2)
    l3 = np.einsum('ijk,j,k->i', C, alpha, alpha) / (32 * np.pi**3)
    yukawa_pull = (yt**2 / (8 * np.pi**2)) * c_t
    
    d_yt = (yt / (16 * np.pi**2)) * (4.5 * yt**2 - 8.0*g_sq[2] - 2.25*g_sq[1] - 0.85*g_sq[0])
    return np.concatenate((-(l1 + l2 + l3 - yukawa_pull), [d_yt]))

def event_landau_pole(t, y, b, B, C):
    """Terminates if any gauge coupling becomes too strong (a_inv < 1.0)."""
    return np.min(y[:3]) - 1.0
event_landau_pole.terminal = True

def event_yukawa_blowup(t, y, b, B, C):
    """Terminates if top Yukawa exits perturbative regime (yt > 3.5)."""
    return 3.5 - y[3]
event_yukawa_blowup.terminal = True

EVENTS = [event_landau_pole, event_yukawa_blowup]
TOLS = {'rtol': 1e-10, 'atol': 1e-12}

# ==============================================================================
# 4. SINGLE UNIVERSE EVALUATOR
# ==============================================================================
def evaluate_universe(s2w_guess, mz, a_em_inv_mz, v_ew, m_top_pole, as_mz, m_warden, m_melt):
    as_pi = as_mz / np.pi
    mt_msbar = m_top_pole * (1.0 - (4.0/3.0)*as_pi - 1.0414*(as_pi**2) - 3.3714*(as_pi**3))
    yt_start = np.sqrt(2) * mt_msbar / v_ew

    phase_width = np.log(m_melt / m_warden)
    Theta_S = -(warden_weights / (2 * np.pi)) * phase_width
    Theta_Theory = Theta_S + Theta_D 

    y_mz = np.array([
        (3/5)*a_em_inv_mz*(1-s2w_guess) + delta_sm_mz[0], 
        a_em_inv_mz*s2w_guess + delta_sm_mz[1], 
        1/as_mz + delta_sm_mz[2], 
        yt_start
    ])
    
    # Phase 1: SM
    s1 = solve_ivp(rge_precision, [np.log(mz), np.log(m_warden)], y_mz, args=(b_SM, B_SM, C_SM), 
                   method='Radau', events=EVENTS, **TOLS)
    if s1.status == 1: raise ValueError("Phase 1 Landau Pole")
    
    yw = s1.y[:, -1].copy()
    yw[:3] += delta_warden
    
    # Phase 2: Warden Scalar Phase
    s2 = solve_ivp(rge_precision, [np.log(m_warden), np.log(m_melt)], yw, args=(b_W, B_W, C_W), 
                   method='Radau', events=EVENTS, **TOLS)
    if s2.status == 1: raise ValueError("Phase 2 Landau Pole")
    
    ym = s2.y[:, -1].copy()
    ym[:3] += Theta_Theory
    
    # Root Finder for Intersection
    def cross(ln_scale):
        s3 = solve_ivp(rge_precision, [np.log(m_melt), ln_scale], ym, args=(b_W, B_W, C_W), 
                       method='Radau', events=EVENTS, **TOLS)
        if s3.status == 1: raise ValueError("High Energy Landau Pole")
        return s3.y[0, -1] - s3.y[1, -1]
    
    res = root_scalar(cross, bracket=[np.log(1e15), np.log(1e18)], method='brentq')
    ln_mgut = res.root
    
    s_f = solve_ivp(rge_precision, [np.log(m_melt), ln_mgut], ym, args=(b_W, B_W, C_W), 
                    method='Radau', events=EVENTS, **TOLS)
    gap = s_f.y[1, -1] - s_f.y[2, -1]
    return gap, ln_mgut

# ==============================================================================
# 5. MONTE CARLO EXECUTION
# ==============================================================================
def sample_truncnorm(mean, std, sigmas=3):
    return truncnorm.rvs(-sigmas, sigmas, loc=mean, scale=std)

def main():
    iterations = 500  # Reduced to 500 for a faster initial test run. Increase to 5000 later.
    s2w_results = []
    mgut_results = []
    failed_runs = 0

    print(f"Initiating Strict 3-Loop Monte Carlo Predictor ({iterations} Universes)...")
    np.random.seed(42) 

    for i in range(iterations):
        # Progress indicator
        if i % 50 == 0 and i > 0:
            sys.stdout.write(f"\rRunning... {i}/{iterations} complete.")
            sys.stdout.flush()

        _mz = sample_truncnorm(*mz_input)
        _a_em = sample_truncnorm(*a_em_inv_input)
        _v_ew = sample_truncnorm(*v_ew_input)
        _m_top = sample_truncnorm(*m_top_input)
        _as_mz = sample_truncnorm(*as_mz_input)
        _m_w = sample_truncnorm(*m_warden_input)
        _m_m = sample_truncnorm(*m_melt_input)

        try:
            def objective(x):
                # Returns the gap between alpha_2 and alpha_3 at the MGUT crossing point
                return evaluate_universe(x, _mz, _a_em, _v_ew, _m_top, _as_mz, _m_w, _m_m)[0]
            
            opt_s2w = root_scalar(objective, bracket=[0.225, 0.238], method='brentq').root
            _, opt_ln_mgut = evaluate_universe(opt_s2w, _mz, _a_em, _v_ew, _m_top, _as_mz, _m_w, _m_m)
            
            s2w_results.append(opt_s2w)
            mgut_results.append(np.exp(opt_ln_mgut))
            
        except Exception as e:
            # Catches everything from Landau poles to root-finder bracket failures
            failed_runs += 1 

    # ==============================================================================
    # 6. RESULTS
    # ==============================================================================
    print("\n\n================================================================")
    print("     U(4) STRICT 3-LOOP PREDICTIONS WITH EVENT TERMINATION")
    print("================================================================")
    
    if iterations == 0:
        return

    success_rate = ((iterations - failed_runs) / iterations) * 100
    print(f"Convergence Rate     : {success_rate:.1f}% ({iterations - failed_runs}/{iterations} runs unified)")

    if len(s2w_results) > 0:
        print(f"OUTPUT sin^2 theta_W : {np.mean(s2w_results):.5f} ± {np.std(s2w_results):.5f}")
        print(f"OUTPUT M_GUT Scale   : {np.mean(mgut_results):.3e} ± {np.std(mgut_results):.3e} GeV")
    else:
        print("ERROR: No runs successfully converged. The phase space may be too strict.")
    print("================================================================")

if __name__ == "__main__":
    main()
