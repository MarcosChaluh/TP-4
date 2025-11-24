"""
Step-by-step plan (Python section)
1. Load Excel data (user sets path/sheet) and slice sample period.
2. Build level series: Y, C, I≈Absorption−Consumption, K, H (hours*persons engaged), productivity (Y/H), real interest rate.
3. Log-transform (except interest rate) and detrend via HP filter (lambda adjustable).
4. Compute moments on cyclical components: std devs and cross-correlations with GDP at lags/leads (-1,0,+1).
5. Calibrate RBC parameters from data (alpha from labor share if provided, delta, beta, TFP AR(1)).
6. Solve a log-linear RBC policy system for c and h as linear functions of states (k,z); derive k,y,i recursively.
7. Simulate T periods with burn-in for multiple seeds; compute moments for each sim and average across runs.
8. Assemble two tables (data vs model) and print/export.
Self-check: all moments computed on detrended series; correlations use aligned samples; simulation length matches data sample (after burn-in removal).
"""
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Tuple, Dict, List
from statsmodels.tsa.filters.hp_filter import hpfilter
from statsmodels.tsa.ar_model import AutoReg

# -----------------------------
# Data handling utilities
# -----------------------------

def load_argentina_data(
    path: str,
    sheet_name: str = 0,
    date_col: str = None,
    start_date: str = None,
    end_date: str = None,
    freq: str = "A"  # adjust if quarterly ("Q")
) -> pd.DataFrame:
    """Load Excel data, set a date index, and subset the sample period."""
    df = pd.read_excel(path, sheet_name=sheet_name)
    if date_col and date_col in df.columns:
        df[date_col] = pd.to_datetime(df[date_col])
        df = df.set_index(date_col).sort_index()
    else:
        # create a simple period index if not provided
        df.index = pd.period_range(start=0, periods=len(df), freq=freq)
    if start_date:
        df = df[df.index >= pd.Period(start_date, freq=freq)]
    if end_date:
        df = df[df.index <= pd.Period(end_date, freq=freq)]
    return df

@dataclass
class MacroSeries:
    Y: pd.Series
    C: pd.Series
    I: pd.Series
    K: pd.Series
    H: pd.Series
    prod: pd.Series
    r: pd.Series


def build_macro_series(df: pd.DataFrame, lambda_hp: float = 1600.0) -> Tuple[MacroSeries, Dict[str, pd.Series]]:
    """Construct macro series and their cyclical (HP-filtered) components."""
    # Extract raw levels
    Y = df["National Accounts-Based Variables, GDP at National Prices, Constant Prices"].astype(float)
    C = df["National Accounts-Based Variables, Real Consumption at National Prices, Constant Prices"].astype(float)
    absorption = df["National Accounts-Based Variables, Real Domestic Absorption at National Prices, Constant Prices"].astype(float)
    I = absorption - C  # approximation, documented
    K = df["National Accounts-Based Variables, Capital Stock at Constant 2017 National Prices, Constant Prices"].astype(float)
    hours = df["Real GDP, Employment & Population Levels, Average Annual Hours Worked by Persons Engaged"].astype(float)
    persons = df["Real GDP, Employment & Population Levels, Number of Persons Engaged"].astype(float)
    H = hours * persons
    prod = Y / H
    r = df["National Accounts-Based Variables, Real Internal Rate of Return"].astype(float)

    # Log-transform where standard
    log_vars = {"Y": np.log(Y), "C": np.log(C), "I": np.log(I), "K": np.log(K), "H": np.log(H), "prod": np.log(prod)}

    # HP filter cyclical component
    cyc = {}
    for name, series in log_vars.items():
        cycle, _ = hpfilter(series, lamb=lambda_hp)
        cyc[name] = cycle
    # Real interest rate detrending (demean only)
    cyc["r"] = r - r.mean()

    macro = MacroSeries(Y=Y, C=C, I=I, K=K, H=H, prod=prod, r=r)
    return macro, cyc

# -----------------------------
# Moment computation
# -----------------------------

def std_corr_table(cyc: Dict[str, pd.Series], var_order: List[str]) -> pd.DataFrame:
    """Compute std dev and cross-correlations with GDP at leads/lags -1,0,+1."""
    results = []
    y = cyc["Y"]
    for var in var_order:
        x = cyc[var]
        sigma = x.std(ddof=1)
        corr_lag = x.corr(y.shift(1))  # corr(x_t, y_{t-1})
        corr_contemp = x.corr(y)
        corr_lead = x.corr(y.shift(-1))  # corr(x_t, y_{t+1})
        results.append({
            "Variable": var,
            "StdDev": sigma,
            "Corr_Y_t-1": corr_lag,
            "Corr_Y_t": corr_contemp,
            "Corr_Y_t+1": corr_lead,
        })
    table = pd.DataFrame(results).set_index("Variable")
    return table

# -----------------------------
# Calibration helpers
# -----------------------------

@dataclass
class Calibration:
    alpha: float
    delta: float
    beta: float
    rho: float
    sigma_eps: float
    theta: float = 2.0  # labor disutility scale
    phi: float = 1.0    # Frisch elasticity parameter (inverse)


def calibrate_parameters(df: pd.DataFrame, cyc: Dict[str, pd.Series], labor_share_col: str = None) -> Calibration:
    """Calibrate alpha, delta, beta, and TFP AR(1) parameters from data."""
    # Capital share
    if labor_share_col and labor_share_col in df.columns:
        labor_share = df[labor_share_col].dropna()
        alpha = 1.0 - labor_share.mean()
    else:
        alpha = 0.33  # fallback

    delta = df["National Accounts-Based Variables, Average Depreciation Rate of the Capital Stock"].astype(float).mean()
    r_mean = df["National Accounts-Based Variables, Real Internal Rate of Return"].astype(float).mean()
    beta = 1.0 / (1.0 + r_mean)

    tfp = np.log(df["National Accounts-Based Variables, Total Factor Productivity at Constant National Prices (2017=1), Constant Prices"].astype(float)).dropna()
    ar1 = AutoReg(tfp.diff().dropna(), lags=1, old_names=False).fit()
    rho = 1 + ar1.params[1]  # because we estimated on differences
    sigma_eps = ar1.resid.std(ddof=1)

    return Calibration(alpha=alpha, delta=delta, beta=beta, rho=rho, sigma_eps=sigma_eps)

# -----------------------------
# RBC linear policy (manual log-linearization)
# -----------------------------

def steady_state(cal: Calibration) -> Dict[str, float]:
    """Compute deterministic steady state for the RBC model."""
    alpha, delta, beta, theta, phi = cal.alpha, cal.delta, cal.beta, cal.theta, cal.phi
    r_k = (1.0 / beta) - 1.0 + delta
    k_over_h = (alpha / r_k) ** (1.0 / (1.0 - alpha))
    h_ss = ((1 - alpha) / theta) ** (1.0 / (phi + alpha))
    y_ss = k_over_h ** alpha * h_ss ** (1 - alpha)
    k_ss = k_over_h * h_ss
    i_ss = delta * k_ss
    c_ss = y_ss - i_ss
    return {"r_k": r_k, "k_ss": k_ss, "h_ss": h_ss, "y_ss": y_ss, "c_ss": c_ss, "i_ss": i_ss}


def solve_linear_policies(cal: Calibration, ss: Dict[str, float]) -> Dict[str, float]:
    """
    Solve for linear policy coefficients for c_hat and h_hat as functions of states (k_hat, z).
    Uses log-linearized resource constraint, capital accumulation, labor FOC, and Euler equation.
    """
    alpha, delta, beta, theta, phi = cal.alpha, cal.delta, cal.beta, cal.theta, cal.phi
    c_share = ss["c_ss"] / ss["y_ss"]
    i_share = ss["i_ss"] / ss["y_ss"]
    r_k = ss["r_k"]
    phi_euler = (alpha * ss["y_ss"] / ss["k_ss"]) / (alpha * ss["y_ss"] / ss["k_ss"] + 1 - delta)

    # Unknown coefficients: c_hat = a1*k_hat + a2*z; derive h_hat from labor FOC
    a1, a2 = sym_solve_coefficients(alpha, delta, phi, c_share, i_share, phi_euler, cal.rho)
    return {"c_k": a1, "c_z": a2, "c_share": c_share, "i_share": i_share, "phi_euler": phi_euler}


def sym_solve_coefficients(alpha, delta, phi, c_share, i_share, phi_euler, rho):
    """Solve for policy coefficients a1, a2 using linear algebra on Euler equation."""
    # Helper functions of a1,a2 -> implied coefficients in Euler equation
    def implied_system(a1, a2):
        def h_coeffs(k_coeff, z_coeff):
            return (
                (alpha - k_coeff) / (phi + alpha),
                (1.0 - z_coeff) / (phi + alpha)
            )
        hk, hz = h_coeffs(a1, a2)
        yk = alpha + (1 - alpha) * hk
        yz = 1.0 + (1 - alpha) * hz
        # investment and capital transitions
        ik = (yk - c_share * a1) / i_share
        iz = (yz - c_share * a2) / i_share
        kk = (1 - delta) + delta * ik
        kz = delta * iz
        # expectations next period
        c_next_k = a1 * kk
        c_next_z = a1 * kz + a2 * rho
        h_next_k = (alpha - c_next_k) / (phi + alpha)
        h_next_z = (rho - c_next_z) / (phi + alpha)
        y_next_k = alpha * kk + (1 - alpha) * h_next_k
        y_next_z = rho + (1 - alpha) * h_next_z
        euler_k = c_next_k - phi_euler * (y_next_k - kk)
        euler_z = c_next_z - phi_euler * (y_next_z - kz)
        return euler_k, euler_z

    # Solve linear system analytically via numpy
    # We use small numerical search; dimension 2 allows closed form via linearization around guess
    A = np.zeros((2, 2))
    b = np.zeros(2)
    # Finite difference to get Jacobian at zero and assume linearity (good for log-linear model)
    base = np.array(implied_system(0.0, 0.0))
    eps = 1e-4
    for i, (dk, dz) in enumerate([(eps, 0.0), (0.0, eps)]):
        val = np.array(implied_system(dk, dz))
        A[:, i] = (val - base) / eps
    b = -np.array(base)
    sol = np.linalg.solve(A, b)
    return sol[0], sol[1]

# -----------------------------
# Simulation
# -----------------------------

def simulate_rbc(cal: Calibration, ss: Dict[str, float], pol: Dict[str, float], T: int, burnin: int = 50, seed: int = 0) -> pd.DataFrame:
    """Simulate the log-linear RBC for T periods (post burn-in)."""
    alpha, delta, rho, phi = cal.alpha, cal.delta, cal.rho, cal.phi
    c_share, i_share = pol["c_share"], pol["i_share"]
    a1, a2 = pol["c_k"], pol["c_z"]
    # Initialize states in log deviations
    k_hat = 0.0
    z = 0.0
    out = []
    rng = np.random.default_rng(seed)
    for t in range(T + burnin):
        eps = rng.normal(scale=cal.sigma_eps)
        z = rho * z + eps
        # policies
        c_hat = a1 * k_hat + a2 * z
        h_hat = (z + alpha * k_hat - c_hat) / (cal.phi + alpha)
        y_hat = z + alpha * k_hat + (1 - alpha) * h_hat
        i_hat = (y_hat - c_share * c_hat) / i_share
        k_hat_next = (1 - delta) * k_hat + delta * i_hat
        r_hat = alpha * (y_hat - k_hat)  # approx real interest log deviation

        # Convert to levels
        c = ss["c_ss"] * np.exp(c_hat)
        h = ss["h_ss"] * np.exp(h_hat)
        y = ss["y_ss"] * np.exp(y_hat)
        i = ss["i_ss"] * np.exp(i_hat)
        k = ss["k_ss"] * np.exp(k_hat)
        r = ss["r_k"] * np.exp(r_hat)
        prod = y / h

        if t >= burnin:
            out.append({"Y": y, "C": c, "I": i, "K": k, "H": h, "prod": prod, "r": r})
        k_hat = k_hat_next
    return pd.DataFrame(out)


def averaged_model_moments(cal: Calibration, data_len: int, n_sim: int = 20, lambda_hp: float = 1600.0) -> pd.DataFrame:
    ss = steady_state(cal)
    pol = solve_linear_policies(cal, ss)
    tables = []
    for seed in range(n_sim):
        sim = simulate_rbc(cal, ss, pol, T=data_len, burnin=50, seed=seed)
        _, cyc = build_macro_series(sim, lambda_hp=lambda_hp)
        table = std_corr_table(cyc, ["Y", "C", "I", "K", "H", "prod"])
        tables.append(table)
    avg_table = pd.concat(tables).groupby(level=0).mean()
    return avg_table

# -----------------------------
# Main runnable example
# -----------------------------

def main():
    # User-specified inputs
    excel_path = "Arg GDP data.xlsx"  # update as needed
    sheet = 0
    start_date = None
    end_date = None
    lambda_hp = 1600.0

    df = load_argentina_data(excel_path, sheet_name=sheet, start_date=start_date, end_date=end_date, freq="A")
    macro, cyc = build_macro_series(df, lambda_hp=lambda_hp)

    var_order = ["Y", "C", "I", "K", "H", "prod"]
    data_table = std_corr_table(cyc, var_order)
    print("\nTable 14.1 style statistics (Data, Argentina):")
    print(data_table)

    cal = calibrate_parameters(df, cyc)
    model_table = averaged_model_moments(cal, data_len=len(df), n_sim=20, lambda_hp=lambda_hp)
    print("\nTable 14.2 style statistics (Model, averaged across simulations):")
    print(model_table)

    # Optional: export to LaTeX
    # print(data_table.to_latex(float_format=lambda x: f"{x:.3f}"))
    # print(model_table.to_latex(float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()

# Summary of outputs:
# - data_table corresponds to Table 14.1 (empirical moments).
# - model_table corresponds to Table 14.2 (RBC-simulated moments).
