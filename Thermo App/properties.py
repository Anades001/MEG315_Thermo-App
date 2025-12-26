import CoolProp.CoolProp as CP

# -------------------------------
# BASIC PROPERTIES
# -------------------------------
def specific_volume(fluid, T, P):
    """Specific volume from T and P"""
    return 1 / CP.PropsSI("D", "T", T, "P", P, fluid)


def internal_energy(fluid, T, P):
    return CP.PropsSI("U", "T", T, "P", P, fluid)

def get_initial_entropy(fluid, P1, x):
    # 'Q' is the CoolProp key for Quality/Dryness Factor
    s1 = CP.PropsSI('S', 'P', P1, 'Q', x, fluid)
    return s1

def entropy_from_v(fluid, P, v):
    """
    Calculates the specific entropy (s) given the fluid, pressure (P), 
    and specific volume (v).
    
    Args:
        fluid (str): The name of the working fluid (e.g., 'Water', 'Air').
        P (float): Pressure in Pascals (Pa).
        v (float): Specific volume in cubic meters per kilogram (m³/kg).
        
    Returns:
        float: Specific entropy in Joules per kilogram Kelvin (J/kg/K).
    """
    
    # Specific volume (v) is the inverse of density (D).
    D = 1.0 / v
    
    # Use CoolProp's PropsSI to look up the specific entropy ('S') 
    # based on Pressure ('P') and Density ('D').
    s = CP.PropsSI('S', 'P', P, 'D', D, fluid)
    
    return s


def entropy_TP(fluid, T, P):
    """Entropy from T and P (single-phase only)"""
    return CP.PropsSI("S", "T", T, "P", P, fluid)


def entropy_from_Pv(fluid, P, v):
    """Entropy from pressure and specific volume
       Safe for steam in two-phase region
    """
    return CP.PropsSI("S", "P", P, "D", 1/v, fluid)


def temperature_from_PS(fluid, P, s):
    return CP.PropsSI("T", "P", P, "S", s, fluid)
