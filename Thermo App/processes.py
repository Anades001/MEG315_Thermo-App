import numpy as np
import CoolProp.CoolProp as CP

from properties import (
    specific_volume,
    entropy_from_Pv,
    entropy_from_v,
    temperature_from_PS,
)

# --- processes.py (Modified constant_volume function) ---
import CoolProp.CoolProp as CP
import numpy as np

# --- In processes.py ---
def constant_volume(fluid, T1_k, P1_pa, P2_pa, V_total, m, steps=50):
    v_spec = V_total / m
    dens = 1.0 / v_spec
    
    # Initial states
    u1 = CP.PropsSI("U", "T", T1_k, "D", dens, fluid)
    s1 = CP.PropsSI("S", "T", T1_k, "D", dens, fluid)
    
    # Final states
    T2_k = CP.PropsSI("T", "P", P2_pa, "D", dens, fluid)
    u2 = CP.PropsSI("U", "P", P2_pa, "D", dens, fluid)
    s2 = CP.PropsSI("S", "P", P2_pa, "D", dens, fluid)
    
    # Numerical results
    heat = m * (u2 - u1) 
    delta_s = m * (s2 - s1)

    # Plotting arrays
    Ps = np.linspace(P1_pa, P2_pa, steps)
    Ts, Ss, Vs = [], [], []
    for P in Ps:
        Ts.append(CP.PropsSI("T", "P", P, "D", dens, fluid))
        Ss.append(CP.PropsSI("S", "P", P, "D", dens, fluid))
        Vs.append(V_total)

    # MUST RETURN 7 ITEMS
    return Ts, Ps, Vs, Ss, T2_k, heat, delta_s


# ------------------------------------------------
# 2. CONSTANT PRESSURE (ISOBARIC)
# ------------------------------------------------
import CoolProp.CoolProp as CP
import numpy as np

def constant_pressure(fluid, T1_k, P_pa, m, steps=50):
    # 1. State 1 Properties
    # Get initial density (D) to find specific volume (v1) and enthalpy (h1)
    d1 = CP.PropsSI("D", "T", T1_k, "P", P_pa, fluid)
    v1 = 1.0 / d1
    h1 = CP.PropsSI("H", "T", T1_k, "P", P_pa, fluid)

    # 2. State 2 Properties
    # The question states the volume is doubled
    v2 = v1 * 2
    d2 = 1.0 / v2
    
    # Find final temperature (T2) and enthalpy (h2) at P and v2
    T2_k = CP.PropsSI("T", "P", P_pa, "D", d2, fluid)
    h2 = CP.PropsSI("H", "P", P_pa, "D", d2, fluid)

    # 3. Numerical Solutions
    # Work W = P * (V2 - V1) = P * m * (v2 - v1)
    work = P_pa * m * (v2 - v1)
    
    # Heat Q = m * (h2 - h1) for constant pressure
    heat = m * (h2 - h1)

    # 4. Generate arrays for Plotting
    Ts_plot = np.linspace(T1_k, T2_k, steps)
    Vs, Ss, Ps = [], [], []

    for T in Ts_plot:
        d_val = CP.PropsSI("D", "T", T, "P", P_pa, fluid)
        s_val = CP.PropsSI("S", "T", T, "P", P_pa, fluid)
        Vs.append((1.0 / d_val) * m)
        Ss.append(s_val)
        Ps.append(P_pa)

    # Return: Plot data + Numerical answers
    return Ts_plot, Ps, Vs, Ss, T2_k, work, heat


# Polytropic processes

import numpy as np

def polytropic(fluid, T1_k, P1_pa, P2_pa, V1_input, m, n, x1):
    # --- 1. Constants and Initial State ---
    if fluid == "Air":
        R = 287
        cv = 718
    elif fluid == "Methane":
        R = 518
        cv = 1735
    else: # Steam (Simplified for the assignment logic)
        R = 461.5
        cv = 1410 

    # If V1 isn't provided, calculate it from Ideal Gas Law (or dryness for steam)
    # For steam: V1 = x * m * R * T1 / P1
    V1 = (x1 * m * R * T1_k) / P1_pa if V1_input == 0 else V1_input

    # --- 2. Calculate Final State ---
    # T2 = T1 * (P2/P1)^((n-1)/n)
    T2_k = T1_k * (P2_pa / P1_pa)**((n - 1) / n)
    # V2 = V1 * (P1/P2)^(1/n)
    V2 = V1 * (P1_pa / P2_pa)**(1/n)

    # --- 3. Work, Internal Energy, and Heat ---
    # Work W = (P1V1 - P2V2) / (n-1)
    work = (P1_pa * V1 - P2_pa * V2) / (n - 1)
    
    # Change in Internal Energy dU = m * cv * (T2 - T1)
    du = m * cv * (T2_k - T1_k)
    
    # Heat Transfer Q = W + dU
    heat = work + du
    

    # --- 4. Plotting Data ---
    # Generate arrays for smooth curves on the graphs
    P_coords = np.linspace(P1_pa, P2_pa, 50)
    V_coords = V1 * (P1_pa / P_coords)**(1/n)
    T_coords = T1_k * (P_coords / P1_pa)**((n - 1) / n)
    # Simplified entropy for plotting
    S_coords = np.linspace(0, (heat/T1_k), 50) 

    return T_coords, P_coords, V_coords, S_coords, work, heat
# ------------------------------------------------

# Isentropic
# --- Inside processes.py ---
import CoolProp.CoolProp as CP
import numpy as np

def isentropic(fluid, T1_k, P1_pa, P2_pa, V1, m, steps=50):
    # 1. Get initial properties (State 1)
    # s1 is constant for isentropic, u1 is used for work
    s1 = CP.PropsSI("S", "T", T1_k, "P", P1_pa, fluid)
    u1 = CP.PropsSI("U", "T", T1_k, "P", P1_pa, fluid)

    # 2. Get final properties (State 2: P2 and s2 = s1)
    T2_k = CP.PropsSI("T", "P", P2_pa, "S", s1, fluid)
    u2 = CP.PropsSI("U", "P", P2_pa, "S", s1, fluid)
    
    # 3. Calculate Numerical Solution for the assignment
    # Work (W) = m * (u1 - u2) for an adiabatic process (Q=0)
    work = m * (u1 - u2)

    # 4. Generate arrays for the graphs
    Ps = np.linspace(P1_pa, P2_pa, steps)
    Ts, Vs, Ss = [], [], []

    for P in Ps:
        T_val = CP.PropsSI("T", "P", P, "S", s1, fluid)
        # Specific volume v = 1 / density
        v_val = 1.0 / CP.PropsSI("D", "P", P, "S", s1, fluid)
        
        Ts.append(T_val)
        Vs.append(v_val * m)
        Ss.append(s1)

    # Return everything to app.py
    return Ts, Ps, Vs, Ss, work, T2_k
# ------------------------------------------------
# 5. ISOTHERMAL PROCESS
import CoolProp.CoolProp as CP
import numpy as np

def isothermal(fluid, T_k, P1_pa, v2_target, m, steps=50):
    # 1. State 1 (Initial)
    # Get initial density to calculate u1 and s1
    d1 = CP.PropsSI("D", "T", T_k, "P", P1_pa, fluid)
    u1 = CP.PropsSI("U", "T", T_k, "D", d1, fluid)
    s1 = CP.PropsSI("S", "T", T_k, "D", d1, fluid)
    v1 = 1.0 / d1

    # 2. State 2 (Final - using v2=0.28 converted to density)
    d2 = 1.0 / v2_target # Density is the reciprocal of specific volume
    
    # Use "D" instead of "V" to avoid the ValueError
    u2 = CP.PropsSI("U", "T", T_k, "D", d2, fluid)
    s2 = CP.PropsSI("S", "T", T_k, "D", d2, fluid)
    P2_pa = CP.PropsSI("P", "T", T_k, "D", d2, fluid)

    # 3. Numerical Calculations
    delta_u = m * (u2 - u1)      # Change in Internal Energy
    delta_s = m * (s2 - s1)      # Change in Entropy
    heat = T_k * delta_s         # Q = T * ΔS
    work = heat - delta_u        # W = Q - ΔU

    # 4. Arrays for Plotting
    # Generate a range of densities from d1 to d2
    ds_plot = np.linspace(d1, d2, steps)
    Ps, Vs, Ss, Ts = [], [], [], []

    for d in ds_plot:
        Ps.append(CP.PropsSI("P", "T", T_k, "D", d, fluid))
        Vs.append((1.0 / d) * m)
        Ss.append(CP.PropsSI("S", "T", T_k, "D", d, fluid))
        Ts.append(T_k)

    return Ts, Ps, Vs, Ss, delta_u, delta_s, heat, work


#Throttling process
def throttling_process(fluid, P1_bar, T1_C, P2_bar, mass_flow_g_s, diameter_cm):
    # Convert units
    P1 = P1_bar * 1e5
    P2 = P2_bar * 1e5
    T1 = T1_C + 273.15
    m_dot = mass_flow_g_s / 1000 # kg/s
    area = np.pi * (diameter_cm / 200)**2

    # 1. State 1
    h1 = CP.PropsSI('H', 'P', P1, 'T', T1, fluid)
    s1 = CP.PropsSI('S', 'P', P1, 'T', T1, fluid)
    v1 = 1.0 / CP.PropsSI('D', 'P', P1, 'T', T1, fluid)

    # 2. State 2 (Isenthalpic: h2 = h1)
    # Throttling is a constant enthalpy process
    T2 = CP.PropsSI('T', 'P', P2, 'H', h1, fluid)
    s2 = CP.PropsSI('S', 'P', P2, 'H', h1, fluid)
    v2 = 1.0 / CP.PropsSI('D', 'P', P2, 'H', h1, fluid)

    # 3. Calculate Velocity (V = m_dot * v / Area)
    velocity2 = (m_dot * v2) / area
    ds = s2 - s1 # Entropy change

    # Plotting arrays
    Ps = np.linspace(P1, P2, 50)
    Hs = np.full_like(Ps, h1)
    Ts = [CP.PropsSI('T', 'P', p, 'H', h1, fluid) for p in Ps]
    Ss = [CP.PropsSI('S', 'P', p, 'H', h1, fluid) for p in Ps]
    Vs = [(1.0 / CP.PropsSI('D', 'P', p, 'H', h1, fluid)) for p in Ps]

    return Ts, Ps, Vs, Ss, velocity2, ds


#Steam vessel filling
def steam_vessel_filling(V_vessel, P1_bar, x1, P_main_bar, T_main_C):
    # Units conversion
    P1 = P1_bar * 1e5
    P_main = P_main_bar * 1e5
    T_main = T_main_C + 273.15
    
    # 1. Initial State in Vessel
    # Use Quality (x) and Pressure (P) to find initial density and internal energy
    d1 = CP.PropsSI('D', 'P', P1, 'Q', x1, 'Water')
    u1 = CP.PropsSI('U', 'P', P1, 'Q', x1, 'Water')
    m1 = d1 * V_vessel
    
    # 2. State of Steam in Main (Inlet enthalpy h_in)
    h_in = CP.PropsSI('H', 'P', P_main, 'T', T_main, 'Water')
    
    # 3. Final State (Energy Balance: m2*u2 - m1*u1 = (m2-m1)*h_in)
    # This requires an iterative solve or property lookup where u2 + (m1/m2)(h_in - u1) = h_in
    # For a simplified estimate at the main pressure:
    T2 = CP.PropsSI('T', 'P', P_main, 'H', h_in, 'Water') # Approximation for filling
    d2 = CP.PropsSI('D', 'P', P_main, 'T', T2, 'Water')
    m2 = d2 * V_vessel
    
    return m2, T2, (m2 - m1)


#Turbine discharge
def turbine_discharge(V_tank, P1_bar, T1_C, P2_bar):
    # Units
    P1 = P1_bar * 1e5
    P2 = P2_bar * 1e5
    T1 = T1_C + 273.15
    
    # 1. Initial State
    d1 = CP.PropsSI('D', 'P', P1, 'T', T1, 'Air')
    u1 = CP.PropsSI('U', 'P', P1, 'T', T1, 'Air')
    s1 = CP.PropsSI('S', 'P', P1, 'T', T1, 'Air')
    m_initial = d1 * V_tank
    
    # 2. Final State (Isentropic expansion: s2 = s1)
    u2 = CP.PropsSI('U', 'P', P2, 'S', s1, 'Air')
    T2 = CP.PropsSI('T', 'P', P2, 'S', s1, 'Air')
    
    # 3. Maximum Work = m * (u1 - u2)
    work_max = m_initial * (u1 - u2)
    
    # Plotting arrays
    Ps = np.linspace(P1, P2, 50)
    Ts = [CP.PropsSI('T', 'P', p, 'S', s1, 'Air') for p in Ps]
    Ss = np.full_like(Ps, s1) # Vertical line on T-s diagram
    Vs = [(m_initial / CP.PropsSI('D', 'P', p, 'S', s1, 'Air')) for p in Ps]
    
    return Ts, Ps, Vs, Ss, work_max, m_initial