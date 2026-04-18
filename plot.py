import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

# ==============================================================================
# 1. CENTRAL EXPERIMENTAL INPUTS
# ==============================================================================
mz = 91.1876
a_em_inv_mz = 127.955
v_ew = 246.21965
m_top_pole = 172.57
as_mz = 0.1180

m_warden = 8165.86  
m_melt = 257608.53  

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
# 3. ODE SOLVER
# ==============================================================================
def rge_precision(t, y, b, B, C):
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

# ==============================================================================
# 4. TRAJECTORY GENERATOR
# ==============================================================================
def get_trajectories(s2w_guess):
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
    
    s1 = solve_ivp(rge_precision, [np.log(mz), np.log(m_warden)], y_mz, args=(b_SM, B_SM, C_SM), 
                   method='Radau', rtol=1e-10, atol=1e-12, dense_output=True)
    
    yw = s1.y[:, -1].copy()
    yw[:3] += delta_warden
    
    s2 = solve_ivp(rge_precision, [np.log(m_warden), np.log(m_melt)], yw, args=(b_W, B_W, C_W), 
                   method='Radau', rtol=1e-10, atol=1e-12, dense_output=True)
    
    ym = s2.y[:, -1].copy()
    ym[:3] += Theta_Theory
    
    def cross(ln_scale):
        s3 = solve_ivp(rge_precision, [np.log(m_melt), ln_scale], ym, args=(b_W, B_W, C_W), 
                       method='Radau', rtol=1e-10, atol=1e-12)
        return s3.y[0, -1] - s3.y[1, -1]
    
    res = root_scalar(cross, bracket=[np.log(1e15), np.log(1e18)], method='brentq')
    ln_mgut = res.root
    
    s3 = solve_ivp(rge_precision, [np.log(m_melt), np.log(1e19)], ym, args=(b_W, B_W, C_W), 
                   method='Radau', rtol=1e-10, atol=1e-12, dense_output=True)
                   
    return s1, s2, s3, ln_mgut

def find_central_s2w():
    def objective(x):
        s1, s2, s3, ln_mgut = get_trajectories(x)
        y_gut = s3.sol(ln_mgut)
        return y_gut[1] - y_gut[2]
    return root_scalar(objective, bracket=[0.225, 0.238], method='brentq').root

print("Calculating precision trajectories...")
optimal_s2w = find_central_s2w()
s1, s2, s3, ln_mgut = get_trajectories(optimal_s2w)
gut_scale = np.exp(ln_mgut)
alpha_gut_inv = s3.sol(ln_mgut)[0]

# ==============================================================================
# 5. PLOTTING WITH EXPLICIT GEOMETRIC SPIN-TRANSITION
# ==============================================================================
# Analytical derivation: Delta_GUT = (11/3) * (C_2(U4) / pi) = (11 * 4) / (3 * pi)
SPIN_TRANSITION_SHIFT = 44 / (3 * np.pi) 
native_gut_inv = alpha_gut_inv + SPIN_TRANSITION_SHIFT

print(f"Mathematical MS-bar intersection : {alpha_gut_inv:.2f}")
print(f"Geometric Spin-Transition Shift    : +{SPIN_TRANSITION_SHIFT:.3f} (44/3pi)")
print(f"Native U(4) Geometric Coupling     : {native_gut_inv:.1f}")

t1 = np.linspace(np.log(mz), np.log(m_warden), 100)
t2 = np.linspace(np.log(m_warden), np.log(m_melt), 100)
t3 = np.linspace(np.log(m_melt), np.log(1e19), 300)

y1 = s1.sol(t1)
y2 = s2.sol(t2)
y3 = s3.sol(t3)

mu = np.exp(np.concatenate((t1, t2, t3)))
a1_inv = np.concatenate((y1[0], y2[0], y3[0]))
a2_inv = np.concatenate((y1[1], y2[1], y3[1]))
a3_inv = np.concatenate((y1[2], y2[2], y3[2]))

# --- Generate the Plot ---
plt.figure(figsize=(11, 7), dpi=300)

plt.plot(mu, a1_inv, 'k-', linewidth=1.5, label=r'$\alpha_1^{-1}$ (U(1)$_Y$, Hypercharge)')
plt.plot(mu, a2_inv, 'k-', linewidth=1.5, label=r'$\alpha_2^{-1}$ (SU(2)$_L$, Weak)')
plt.plot(mu, a3_inv, 'k-', linewidth=1.5, label=r'$\alpha_3^{-1}$ (SU(3)$_C$, Strong)')

# Melting Threshold Marker
plt.axvline(x=m_melt, color='black', linestyle=':', linewidth=1.2)
plt.text(m_melt * 1.5, 45, r'Warden Threshold ($M_{eff} \approx 258$ TeV)', fontsize=10, verticalalignment='center')

# Unification Boundary Markers
plt.axvline(x=gut_scale, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

# Plot MS-Bar Intersection (Hollow Circle)
plt.plot(gut_scale, alpha_gut_inv, 'ko', fillstyle='none', markersize=8, markeredgewidth=1.5)

# Plot Native U(4) Vacuum Intersection (Solid Circle)
plt.plot(gut_scale, native_gut_inv, 'ko', markersize=7)

# Draw the Spin-Transition Arrow
plt.annotate('', xy=(gut_scale, native_gut_inv - 0.6), xytext=(gut_scale, alpha_gut_inv + 0.6),
             arrowprops=dict(arrowstyle="->", color='red', lw=2, shrinkA=0, shrinkB=0))

# Annotations for the GUT Boundary
plt.text(gut_scale * 0.8, alpha_gut_inv, f'$\\overline{{MS}}$ Intersection\n$\\alpha_{{eff}}^{{-1}} \\approx {alpha_gut_inv:.1f}$', 
         fontsize=10, horizontalalignment='right', verticalalignment='center')

plt.text(gut_scale * 1.2, (alpha_gut_inv + native_gut_inv)/2, 
         f'Warden Spin-Transition\n$\\Delta_{{GUT}} = \\frac{{44}}{{3\\pi}} \\approx +{SPIN_TRANSITION_SHIFT:.2f}$', 
         color='red', fontsize=10, horizontalalignment='left', verticalalignment='center')

plt.text(gut_scale * 0.8, native_gut_inv, f'Native $U(4)$ Vacuum\n$M_{{GUT}} \\approx {gut_scale/1e16:.2f} \\times 10^{{16}}$ GeV\n$\\alpha_{{native}}^{{-1}} \\equiv {native_gut_inv:.1f}$', 
         fontsize=11, horizontalalignment='right', verticalalignment='center', fontweight='bold')

# Formatting
plt.xscale('log')
plt.xlim(1e2, 1e19)
plt.ylim(15, 65)
plt.xlabel(r'Energy Scale ($\mu$) [GeV]', fontsize=12)
plt.ylabel(r'Inverse Coupling ($\alpha_i^{-1}$)', fontsize=12)
plt.title('Precision Gauge Coupling Unification in the U(4) Grand Unified Theory', fontsize=14)
plt.legend(loc='lower right')
plt.grid(True, which="both", ls="-", alpha=0.2)

plt.tight_layout()
plt.savefig('FigureW_GUT_Unification.pdf', format='pdf', bbox_inches='tight')
print("Saved as FigureW_GUT_Unification.pdf")
