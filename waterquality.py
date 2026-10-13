"""
River water quality simulation: Streeter-Phelps dissolved oxygen (DO) sag.

Setup:  python -m pip install streamlit numpy pandas
Run:    python -m streamlit run .\\water_quality_sim.py
"""
import numpy as np
import pandas as pd
import streamlit as st

st.set_page_config(page_title="River Water Quality Simulation", layout="wide")
st.title("River water quality simulation")
st.caption("Streeter-Phelps model: how a pollution discharge affects BOD and dissolved oxygen downstream.")

# ---------------- Inputs ----------------
sb = st.sidebar
sb.header("River")
Qr = sb.slider("River flow (m³/s)", 0.5, 100.0, 10.0)
u = sb.slider("Flow velocity (m/s)", 0.05, 2.0, 0.30)
H = sb.slider("Average depth (m)", 0.3, 10.0, 2.0)
T = sb.slider("Water temperature (°C)", 0.0, 35.0, 20.0)
BODr = sb.slider("Upstream BOD (mg/L)", 0.0, 20.0, 2.0)
DOr = sb.slider("Upstream DO (mg/L)", 0.0, 14.0, 8.0)

sb.header("Discharge")
Qw = sb.slider("Discharge flow (m³/s)", 0.0, 20.0, 1.0)
BODw = sb.slider("Discharge BOD (mg/L)", 0.0, 400.0, 150.0)
DOw = sb.slider("Discharge DO (mg/L)", 0.0, 10.0, 1.0)

sb.header("Reaction rates (at 20 °C)")
kd20 = sb.slider("BOD decay rate kd (1/day)", 0.05, 1.0, 0.30)
auto_ka = sb.checkbox("Estimate reaeration rate from river (O'Connor-Dobbins)", True)
ka20_manual = sb.slider("Reaeration rate ka (1/day)", 0.05, 5.0, 0.8, disabled=auto_ka)

sb.header("Analysis")
length = sb.slider("Reach length (km)", 10, 300, 100)
standard = sb.slider("DO standard (mg/L)", 2.0, 8.0, 5.0)
runs = sb.slider("Monte Carlo runs", 50, 2000, 500)
unc = sb.slider("Rate uncertainty (± fraction)", 0.0, 0.6, 0.25)
seed = sb.number_input("Random seed", value=42, step=1)


# ---------------- Model ----------------
def do_saturation(temp):
    """Approximate DO saturation (mg/L) in fresh water at sea level."""
    return 14.652 - 0.41022 * temp + 0.007991 * temp**2 - 0.000077774 * temp**3


def profile(x_km, L0, D0, kd, ka, vel):
    """BOD L and oxygen deficit D at distances x_km downstream of the discharge."""
    t = x_km * 1000.0 / (vel * 86400.0)  # travel time in days
    diff = ka - kd
    diff = np.where(np.abs(diff) < 1e-6, 1e-6, diff)
    L = L0 * np.exp(-kd * t)
    D = kd * L0 / diff * (np.exp(-kd * t) - np.exp(-ka * t)) + D0 * np.exp(-ka * t)
    return L, D


Q = Qr + Qw
L0 = (Qr * BODr + Qw * BODw) / Q  # fully mixed BOD at the discharge point
DO0 = (Qr * DOr + Qw * DOw) / Q  # fully mixed DO at the discharge point
DOsat = do_saturation(T)
D0 = max(DOsat - DO0, 0.0)

ka20 = 3.93 * np.sqrt(u) / H**1.5 if auto_ka else ka20_manual
kd = kd20 * 1.047 ** (T - 20)  # temperature-corrected rates
ka = ka20 * 1.024 ** (T - 20)

x = np.linspace(0, length, 300)
L_det, D_det = profile(x, L0, D0, kd, ka, u)
DO_det = np.clip(DOsat - D_det, 0, None)

rng = np.random.default_rng(int(seed))
kd_s = kd * np.exp(rng.normal(0, unc, (runs, 1)))
ka_s = ka * np.exp(rng.normal(0, unc, (runs, 1)))
_, D_s = profile(x, L0, D0, kd_s, ka_s, u)
DO_s = np.clip(DOsat - D_s, 0, None)

# ---------------- Results ----------------
i_min = int(np.argmin(DO_det))
p_viol = float((DO_s.min(axis=1) < standard).mean())

c1, c2, c3, c4 = st.columns(4)
c1.metric("Mixed BOD at discharge", f"{L0:.1f} mg/L")
c2.metric("Lowest DO (best estimate)", f"{DO_det[i_min]:.2f} mg/L", f"at {x[i_min]:.0f} km", delta_color="off")
c3.metric("Chance DO falls below standard", f"{p_viol:.0%}")
c4.metric("DO saturation at this temperature", f"{DOsat:.2f} mg/L")

if DO_det[i_min] < standard:
    st.error(f"Dissolved oxygen drops below {standard:.1f} mg/L near {x[i_min]:.0f} km downstream.")
else:
    st.success(f"Dissolved oxygen stays above {standard:.1f} mg/L in the best-estimate run.")

idx = pd.Index(np.round(x, 2), name="Distance downstream (km)")
do_df = pd.DataFrame({
    "Best estimate": DO_det,
    "10th percentile": np.percentile(DO_s, 10, axis=0),
    "90th percentile": np.percentile(DO_s, 90, axis=0),
    "Standard": np.full_like(x, standard),
}, index=idx)
bod_df = pd.DataFrame({"BOD (mg/L)": L_det}, index=idx)

left, right = st.columns(2)
left.subheader("Dissolved oxygen")
left.line_chart(do_df)
right.subheader("Biochemical oxygen demand (BOD)")
right.line_chart(bod_df)

st.download_button("Download results (CSV)", do_df.join(bod_df).to_csv(), "water_quality_results.csv", "text/csv")

with st.expander("Model and assumptions"):
    st.markdown(f"""
- Steady-state, one-dimensional river with one point-source discharge, fully mixed at the discharge point.
- BOD decays at rate kd and oxygen is replenished from the air at rate ka; both are temperature-corrected
  (current values: kd = {kd:.2f}/day, ka = {ka:.2f}/day).
- Not included: nitrification, algae and photosynthesis, sediment oxygen demand, tributaries, dispersion.
- Monte Carlo runs vary kd and ka randomly (lognormal) to show how uncertain rates affect the result.
- This is a teaching and screening tool. Calibrate against measured data before using it for decisions.
""")