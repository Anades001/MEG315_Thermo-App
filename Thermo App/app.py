import streamlit as st
from processes import *
from plots import *
import pandas as pd

def export_results(data_dict):
    df = pd.DataFrame([data_dict])
    csv = df.to_csv(index=False).encode('utf-8')
    return csv


st.set_page_config(page_title="Non-Flow Thermodynamics App", layout="wide")

st.title("Non-Flow Thermodynamic Process Calculator")


# ------------------------------------------------
# --- app.py Selection Sidebar ---
# Add "Home" as the first option in the list
# --- Updated Sidebar ---
category = st.sidebar.radio("Analysis Category", ["Closed Systems", "Flow Processes"])

if category == "Closed Systems":
    process = st.sidebar.selectbox("Select Process", ["Home", "Constant Volume", "Constant Pressure", "Isothermal", "Isentropic", "Polytropic"])
else:
    process = st.sidebar.selectbox("Select Flow Process", ["Home", "Throttling Process", "Turbine Discharge", "Steam Vessel Filling"])

# --- Logic for the Home Page ---
if process == "Home":
    st.title("Thermodynamics Process Simulator")
    st.write("Welcome to the Thermodynamic Analysis Tool. Use the sidebar to select a process and solve the assignment questions.")
    
    # Overview Table based on the assignment images
    st.header("Project Overview")
    st.markdown("""
    This app solves for **Final Temperature**, **Work Done**, **Heat Transfer**, and **Entropy Change** for:
    * **Steam** (Pure Substance)
    * **Air** (Ideal Gas)
    * **Methane** (Ideal Gas)
    """)

    # Instructions
    with st.expander("How to use"):
        st.write("1. Select a process from the sidebar.")
        st.write("2. Input the values from the Sample Questions (e.g., P, T, m).")
        st.write("3. Click **'Run Analysis'** to see the numerical solutions and graphs.")

# --- Then continue with your elif blocks ---
elif process == "Constant Volume":
    # Your existing Constant Volume code...

    # --- app.py (Modified Constant Volume section) ---
# --- Correct Structure for app.py ---
    st.subheader("Constant Volume Analysis")
    
    # 1. Inputs must be defined first
    fluid_raw = st.selectbox("Select Fluid", ["Steam", "Air", "Methane"], key="cv_fluid")
    V_fixed = st.number_input("Fixed Volume (m³)", value=0.14, key="V_cv")
    P1 = st.number_input("Initial Pressure (bar)", value=10.0, key="P1_cv")
    T1 = st.number_input("Initial Temperature (°C)", value=250.0, key="T1_cv")
    P2 = st.number_input("Final Pressure (bar)", value=3.5, key="P2_cv")
    m = st.number_input("Mass (kg)", value=1.0, key="m_cv")

    # 2. The Button must be at the SAME indentation level as the inputs
    if st.button("Run Constant Volume Analysis", key="btn_cv"):
        # Everything below this only happens AFTER the button is clicked
        fluid = "Water" if fluid_raw == "Steam" else fluid_raw
        
        Ts, Ps, Vs, Ss, T2_final, heat_val, dS_val = constant_volume(
            fluid, T1 + 273.15, P1 * 1e5, P2 * 1e5, V_fixed, m
        )
        
        # Display Results
        st.metric("Final Temp (T₂)", f"{T2_final - 273.15:.2f} °C")

        # Display Numerical Results
        c1, c2, c3 = st.columns(3)
        c1.metric("Final Temp (T₂)", f"{T2_final - 273.15:.2f} °C")
        c2.metric("Heat Transfer (Q)", f"{heat_val / 1000:.2f} kJ")
        c3.metric("Entropy Change (ΔS)", f"{dS_val:.4f} kJ/K")

        # Graphs
        st.plotly_chart(plot_pv(Ps, Vs), use_container_width=True)
        st.plotly_chart(plot_ts(Ts, Ss), use_container_width=True)

# ------------------------------------------------
elif process == "Constant Pressure":
    st.subheader("Constant Pressure Process Analysis")
    
    fluid_raw = st.selectbox("Select Fluid", ["Steam", "Air", "Methane"], key="cp_fluid")
    fluid = "Water" if fluid_raw == "Steam" else fluid_raw

    # Inputs from Sample Question
    m = st.number_input("Mass (kg)", value=0.2, key="m_cp")
    T1 = st.number_input("Initial Temperature T₁ (°C)", value=165.0, key="T1_cp")
    P = st.number_input("Constant Pressure (bar)", value=7.0, key="P_cp")

    if st.button("Run Constant Pressure Analysis", key="btn_cp"):
        # Unpack the 7 values
        Ts, Ps, Vs, Ss, T2_final, work_val, heat_val = constant_pressure(
            fluid, T1 + 273.15, P * 1e5, m
        )

        # Display Numerical Results
        res1, res2, res3 = st.columns(3)
        res1.metric("Final Temp (T₂)", f"{T2_final - 273.15:.2f} °C")
        res2.metric("Work (W)", f"{work_val / 1000:.2f} kJ")
        res3.metric("Heat (Q)", f"{heat_val / 1000:.2f} kJ")

        # Graphs
        st.plotly_chart(plot_pv(Ps, Vs), use_container_width=True)
        st.plotly_chart(plot_ts(Ts, Ss), use_container_width=True)

# ------------------------------------------------
elif process == "Polytropic":
    st.subheader("Polytropic Process Analysis")
    
    # --- INPUT SECTION ---
    col1, col2 = st.columns(2)
    
    with col1:
        fluid = st.selectbox("Select Fluid", ["Steam", "Air", "Methane"])
        # Defining 'm' (mass) as required by the sample question
        m = st.number_input("Mass (kg)", value=0.9)
        T1 = st.number_input("T₁ (°C)", value=250.0)
        P1 = st.number_input("P₁ (bar)", value=15.0)

    with col2:
        P2 = st.number_input("P₂ (bar)", value=1.5)
        # Defining 'n' (Polytropic Index)
        n = st.number_input("Polytropic Index (n)", value=1.25)
        # Defining 'V1' - default to 0 so the function calculates it
        V1 = st.number_input("Initial Volume V₁ (m³)", value=0.0, help="Leave 0 to calculate from mass")
        
        # Defining 'x1' (Only needed for Steam)
        x1 = 1.0 # Default for Air/Methane
        if fluid == "Steam":
            x1 = st.number_input("Initial Dryness Factor x₁", 0.0, 1.0, 0.7)

    # --- EXECUTION SECTION ---
    if st.button("Run Analysis"):
        try:
            # We must pass the variables we just defined above
            T, P, V, s, work, heat = polytropic(
                fluid, T1 + 273.15, P1 * 1e5, P2 * 1e5, V1, m, n, x1
            )
            
            # Displaying the numbers for your assignment
            st.success(f"Calculation Complete for {fluid}")
            res_col1, res_col2 = st.columns(2)
            res_col1.metric("Work (W)", f"{work/1000:.2f} kJ")
            res_col2.metric("Heat (Q)", f"{heat/1000:.2f} kJ")
            
            # Plotting
            st.plotly_chart(plot_ts(T, s), use_container_width=True)
            st.plotly_chart(plot_pv(P, V), use_container_width=True)
            
        except Exception as e:
            st.error(f"Error: {e}. Please ensure all inputs are filled.")


# ------------------------------------------------
elif process == "Isentropic":
    st.subheader("Isentropic (Adiabatic) Process")
    
    # Inputs (Same as Polytropic)
    fluid = st.selectbox("Select Fluid", ["Water", "Air", "Methane"]) 
    # Note: Use "Water" for Steam in CoolProp
    m = st.number_input("Mass (kg)", value=0.9)
    T1 = st.number_input("T₁ (°C)", value=250.0)
    P1 = st.number_input("P₁ (bar)", value=15.0)
    P2 = st.number_input("P₂ (bar)", value=1.5)

    # --- Inside app.py ---

if process == "Isentropic":
    # Use names that are friendly for the UI
    fluid_selection = st.selectbox("Select Fluid", ["Steam", "Air", "Methane"])
    
    # Map UI names to CoolProp internal names
    fluid_map = {
        "Steam": "Water",
        "Air": "Air",
        "Methane": "Methane"
    }
    fluid = fluid_map[fluid_selection] # This converts "Steam" to "Water"

    if st.button("Calculate Isentropic"):
        # Now 'fluid' will be "Water", and CoolProp will be happy!
        Ts, Ps, Vs, Ss, work, T2_final = isentropic(
            fluid, T1 + 273.15, P1 * 1e5, P2 * 1e5, 0, m
        )

        # --- Display Numerical Results ---
        st.header("Assignment Results")
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Final Temperature (T₂)", f"{T2_final - 273.15:.2f} °C")
        with col2:
            st.metric("Work Done (W)", f"{work / 1000:.2f} kJ")

        # --- Plotting ---
        st.plotly_chart(plot_pv(Ps, Vs), use_container_width=True)
        st.plotly_chart(plot_ts(Ts, Ss), use_container_width=True)
# ------------------------------------------------
# --- app.py (Modified Isothermal section) ---
elif process == "Isothermal":
    st.subheader("Isothermal Process Analysis")
    
    # Fluid Selection
    fluid_raw = st.selectbox("Select Fluid", ["Steam", "Air", "Methane"], key="iso_fluid")
    fluid = "Water" if fluid_raw == "Steam" else fluid_raw

    # Inputs from Sample Question
    m = st.number_input("Mass (kg)", value=1.0, key="m_iso")
    T1 = st.number_input("Temperature T (°C)", value=155.5, key="T_iso")
    P1 = st.number_input("Initial Pressure P₁ (bar)", value=1.0, key="P1_iso")
    v2 = st.number_input("Final Specific Volume v₂ (m³/kg)", value=0.28, key="v2_iso")

    if st.button("Run Isothermal Analysis", key="btn_iso"):
        # Unpack ALL 8 values returned by the updated function
        Ts, Ps, Vs, Ss, dU, dS, Q, W = isothermal(
            fluid, T1 + 273.15, P1 * 1e5, v2, m
        )

        # Display results required by the assignment
        col1, col2 = st.columns(2)
        col1.metric("Change in Internal Energy (ΔU)", f"{dU/1000:.2f} kJ")
        col1.metric("Change in Entropy (ΔS)", f"{dS/1000:.4f} kJ/K")
        col2.metric("Heat Transfer (Q)", f"{Q/1000:.2f} kJ")
        col2.metric("Work Transfer (W)", f"{W/1000:.2f} kJ")

        # Graphs
        st.plotly_chart(plot_pv(Ps, Vs), use_container_width=True)
        st.plotly_chart(plot_ts(Ts, Ss), use_container_width=True)

#Throttling 
elif process == "Throttling Process":
    st.subheader("Throttling Process Analysis")
    
    col1, col2 = st.columns(2)
    with col1:
        P1 = st.number_input("Inlet Pressure (bar)", value=7.0, key="thr_p1")
        T1 = st.number_input("Inlet Temp (°C)", value=95.0, key="thr_t1")
        m_dot = st.number_input("Flow Rate (g/s)", value=2.3, key="thr_mdot")
    with col2:
        P2 = st.number_input("Outlet Pressure (bar)", value=3.5, key="thr_p2")
        d = st.number_input("Pipe Diameter (cm)", value=15.0, key="thr_d")

    if st.button("Run Throttling Analysis"):
     Ts, Ps, Vs, Ss, velocity2, ds = throttling_process("Air", P1, T1, P2, m_dot, d)

    st.success(f"Downstream Velocity: {velocity2:.4f} m/s")
    st.metric("Entropy Change (Δs)", f"{ds:.4f} J/kg·K")
    st.info("Verification: h₁ = h₂ (Isenthalpic process)")

    # Store results for export
    results = {
        "Process": "Throttling",
        "Fluid": "Air",
        "Inlet_P_bar": P1,
        "Outlet_P_bar": P2,
        "Entropy_Change": ds,
        "Final_Velocity": velocity2
    }
    
    # Add a download button
    st.download_button(
        label="Export Data for SQL",
        data=export_results(results),
        file_name="throttling_results.csv",
        mime="text/csv",
    )

    st.plotly_chart(plot_pv(Ps, Vs), use_container_width=True)
    st.plotly_chart(plot_ts(Ts, Ss), use_container_width=True)


#Steam vessel filling
elif process == "Steam Vessel Filling":
    st.subheader("Steam Vessel Filling Analysis")
    # Sample: 0.085 m3, 1 bar, 0.8 dryness -> 20 bar, 260 C
    
    col1, col2 = st.columns(2)
    with col1:
        st.write("**Vessel Initial State**")
        v_vol = st.number_input("Vessel Volume (m³)", value=0.085, key="fill_v")
        p1 = st.number_input("Initial Pressure (bar)", value=1.0, key="fill_p1")
        x1 = st.number_input("Dryness Fraction (x)", value=0.8, key="fill_x1")
    with col2:
        st.write("**Steam Main State**")
        p_m = st.number_input("Main Pressure (bar)", value=20.0, key="fill_pm")
        t_m = st.number_input("Main Temp (°C)", value=260.0, key="fill_tm")

    if st.button("Estimate Final State"):
        m2, t2, m_added = steam_vessel_filling(v_vol, p1, x1, p_m, t_m)
        st.success(f"Final Mass in Vessel: {m2:.4f} kg")
        st.metric("Mass Added", f"{m_added:.4f} kg")
        st.metric("Final Temp", f"{t2 - 273.15:.2f} °C")

        # SQL Export dictionary
        fill_data = {
            "Process": "Vessel Filling",
            "Initial_Mass_kg": m2 - m_added,
            "Final_Mass_kg": m2,
            "Mass_Added_kg": m_added,
            "Final_Temp_C": t2 - 273.15
        }
        
        st.download_button(
            label="Download Vessel Data for SQL",
            data=export_results(fill_data),
            file_name="steam_filling_results.csv",
            mime="text/csv"
        )


#Turbine discharge
elif process == "Turbine Discharge":
    st.subheader("Maximum Work from Turbine Discharge")
    
    # Inputs from Question
    V_tank = st.number_input("Tank Volume (m³)", value=0.3, key="tur_v")
    P1 = st.number_input("Initial Tank Pressure (bar)", value=35.0, key="tur_p1")
    T1 = st.number_input("Initial Tank Temperature (°C)", value=40.0, key="tur_t1")
    P2 = st.number_input("Atmospheric Pressure (bar)", value=1.0, key="tur_p2")

    if st.button("Calculate Maximum Work"):
        Ts, Ps, Vs, Ss, w_max, m_tank = turbine_discharge(V_tank, P1, T1, P2)
        
        st.header(f"Results for {m_tank:.3f} kg of Air")
        st.metric("Maximum Work (W_max)", f"{w_max / 1000:.2f} kJ")
        st.info("Note: Maximum work assumes an isentropic process (s₁ = s₂).")

        # SQL Export
        tur_data = {
            "Process": "Turbine Discharge",
            "Initial_P_bar": P1,
            "Final_P_bar": P2,
            "Total_Mass_kg": m_tank,
            "Max_Work_kJ": w_max / 1000
        }
        
        st.download_button(
            label="Download Turbine Data for SQL",
            data=export_results(tur_data),
            file_name="turbine_discharge_results.csv",
            mime="text/csv"
        )

        st.plotly_chart(plot_pv(Ps, Vs), use_container_width=True)
        st.plotly_chart(plot_ts(Ts, Ss), use_container_width=True)


#Rankine Cycle
elif process == "Rankine Cycle Template":
    st.subheader("Rankine Cycle Flow Analysis")
    
    # Visualizing the flow logic
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.write("### 1. Boiler")
        p_high = st.number_input("P_high (bar)", value=20.0)
        
    with col2:
        st.write("### 2. Turbine")
        st.info("Expansion to P_low")
        
    with col3:
        st.write("### 3. Condenser")
        p_low = st.number_input("P_low (bar)", value=1.0)
        
    with col4:
        st.write("### 4. Pump")
        st.info("Work Input")