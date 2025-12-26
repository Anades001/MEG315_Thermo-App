import plotly.graph_objects as go

def plot_ts(T, s):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=s, y=T, mode="lines"))
    fig.update_layout(
        title="T–s Diagram",
        xaxis_title="Entropy (J/kg·K)",
        yaxis_title="Temperature (K)"
    )
    return fig

def plot_pv(P, V):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=V, y=P, mode="lines"))
    fig.update_layout(
        title="P–V Diagram",
        xaxis_title="Volume (m³)",
        yaxis_title="Pressure (Pa)"
    )
    return fig
