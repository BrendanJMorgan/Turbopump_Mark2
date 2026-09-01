"""
Fuel properties via CoolProp.

Currently: pure n-dodecane as a Jet-A / RP-1 stand-in.
Properties depend on (T, P) and are valid through subcritical liquid,
supercritical, and pseudo-critical regimes.
"""

import CoolProp.CoolProp as CP

class n_dodecane:
    def __init__(self, name="n-Dodecane"):
        self.name = name
        self._AS = CP.AbstractState("HEOS", name)
        self.Tc = CP.PropsSI("Tcrit", name)   # K
        self.Pc = CP.PropsSI("Pcrit", name)   # Pa
        self.M  = CP.PropsSI("M", name)       # kg/mol

    def _state(self, T, P):
        # CoolProp's PT solver is ill-conditioned near the critical point in
        # both T and P. Nudge away from the singular region.
        dT = self.Tc - T
        dP = self.Pc - P
        
        # If we're inside a small box around the critical point, push out
        if abs(dT) < 0.5 and abs(dP) < 0.05 * self.Pc:  # 0.5 K and 5% of Pc
            # Move radially outward to whichever side we're nearer
            if abs(dT) > abs(dP) / self.Pc * self.Tc:  # dominate by T
                T = self.Tc - 0.5 if dT > 0 else self.Tc + 0.5
            else:
                P = self.Pc * 0.95 if dP > 0 else self.Pc * 1.05
        self._AS.update(CP.PT_INPUTS, P, T)

    # Bulk properties at (T [K], P [Pa])
    def rho(self, T, P): self._state(T, P); return self._AS.rhomass()       # kg/m^3
    def cp (self, T, P): self._state(T, P); return self._AS.cpmass()        # J/kg-K
    def cv (self, T, P): self._state(T, P); return self._AS.cvmass()        # J/kg-K
    def mu (self, T, P): self._state(T, P); return self._AS.viscosity()     # Pa-s
    def k  (self, T, P): self._state(T, P); return self._AS.conductivity()  # W/m-K
    def h  (self, T, P): self._state(T, P); return self._AS.hmass()         # J/kg
    def s  (self, T, P): self._state(T, P); return self._AS.smass()         # J/kg-K
    def a  (self, T, P): self._state(T, P); return self._AS.speed_sound()   # m/s

    def Pr(self, T, P):
        self._state(T, P)
        return self._AS.cpmass() * self._AS.viscosity() / self._AS.conductivity()

    # Phase boundary (subcritical only)
    def is_supercritical(self, P):
        return P >= self.Pc

    def Tsat(self, P):
        """Saturation T [K] at P [Pa]. NaN if supercritical."""
        if P >= self.Pc:
            return float("nan")
        self._AS.update(CP.PQ_INPUTS, P, 0.0)
        return self._AS.T()

    def Psat(self, T):
        """Saturation P [Pa] at T [K]. NaN if T >= Tc."""
        if T >= self.Tc:
            return float("nan")
        self._AS.update(CP.QT_INPUTS, 0.0, T)
        return self._AS.p()

    def hvap(self, T):
        """Latent heat [J/kg] at T [K]. Zero above Tc."""
        if T >= self.Tc:
            return 0.0
        self._AS.update(CP.QT_INPUTS, 0.0, T)
        h_l = self._AS.hmass()
        self._AS.update(CP.QT_INPUTS, 1.0, T)
        h_v = self._AS.hmass()
        return h_v - h_l