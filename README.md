# Thermodynamic Process Simulator (MEG315 Assignment)

This Streamlit application visualizes and calculates thermodynamic properties for various processes involving Steam, Air, and Methane. It was developed to fulfill the requirements for MEG315 Assignments at the University of Lagos.

## Features
- **Closed Systems:** Constant Volume, Constant Pressure, Isothermal, Isentropic, and Polytropic analysis.
- **Flow Processes:** Throttling, Turbine Discharge (Max Work), and Steam Vessel Filling.
- **Data Visualization:** Real-time generation of P-v and T-s diagrams using Plotly.
- **SQL Integration:** Exportable CSV data formatted for SQL Server ingestion.

## Technical Logic
- **Throttling:** Modeled as an isenthalpic process ($h_1 = h_2$) to determine downstream velocity and entropy generation.
- **Turbine Discharge:** Uses isentropic expansion ($s_1 = s_2$) to calculate the maximum potential work from a pressurized tank.
- **Property Tables:** Uses the `CoolProp` library for high-accuracy fluid property lookups.

## How to Run
1. Install dependencies: `pip install -r requirements.txt`
2. Run the app: `streamlit run app.py`
