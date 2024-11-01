import pandas as pd
import matplotlib.pyplot as plt
from typing import Callable
from scipy.integrate import solve_ivp

from parameters import *
from curves import m_inf


METHOD = 'LSODA'
RTOL = 1e-5
ATOL = 1e-8


#TODO: Add documentation
class VoltageResponse:
    def __init__(
        self,
        Vlg: Callable,
        dVlgdt: Callable,
        dVsdt: Callable,
        t_max: float = 100000.0
    ):
        self.Vlg = Vlg
        self.dVlgdt = dVlgdt
        self.dVsdt = dVsdt
        self.t_max = t_max
        self.solution = self._solve()

    def create_equation(self) -> Callable:
        def equation(t, y):
            Vj, Vm = y

            Jj = gKJ * m_inf(Vm - Vj, Vh, km) * (Vm - Vj - VK)
            Jf = gKF * m_inf(Vm - self.Vlg(t), Vh, km) * (Vm - self.Vlg(t) - VK)

            dVjdt = cm * (1 - r) / (cm * (1 - r) + cJ) * self.dVlgdt(t) - gJ / (cm * (1 - r) + cJ) * Vj + (1 - r) / (cm * (1 - r) + cJ) * (Jj - Jf) + cJ / (cm * (1 - r) + cJ) * self.dVsdt(t)
            dVmdt = r * dVjdt + (1 - r) * self.dVlgdt(t) - 1 / cm * ((1 - r) * Jf + r * Jj)

            return [dVjdt, dVmdt]

        return equation

    def _solve(self):
        equation = self.create_equation()
        return solve_ivp(
            equation,
            [0.0, self.t_max],
            [0.0, Vrest],
            method=METHOD,
            rtol=RTOL,
            atol=ATOL
        )

    def as_df(self) -> pd.DataFrame:
        t_sec = self.solution.t * 10 ** (-3)
        V_LG = self.Vlg(self.solution.t)
        V_J = self.solution.y[0, :]
        V_M = self.solution.y[1, :]
        V_JM = V_M - V_J
        V_FM = V_M - V_LG
        P_JM = m_inf(V_JM, Vh, km)
        P_FM = m_inf(V_FM, Vh, km)
        IK_J = AM * r * gKJ * P_JM * (V_JM - VK)
        IK_F = AM * (1 - r) * gKF * P_FM * (V_FM - VK)
        I_J = AM * r * gJ * V_J

        df = pd.DataFrame(
            {
                't(s)' : t_sec,
                'VLG(mV)' : V_LG,
                'VJ(mV)' : V_J,
                'VM(mV)' : V_M,
                'VFM(mV)' : V_FM,
                'VJM(mV)' : V_JM,
                'PFM' : P_FM,
                'PJM' : P_JM,
                'IKJ(uA)' : IK_J,
                'IKF(uA)' : IK_F,
                'IJ(uA)' : I_J
            }
        )

        return df

    def plot(self):
        t_sec = self.solution.t * 10 ** (-3)
        V_LG = self.Vlg(self.solution.t)
        V_J = self.solution.y[0, :]
        V_M = self.solution.y[1, :]

        plt.figure(figsize=(12, 6))
        ax = plt.axes((0.1, 0.1, 0.8, 0.8))

        ax.plot(t_sec, V_LG, color='red', zorder = 1)
        p_VJ = ax.plot(t_sec, V_J, color='blue', marker='o', markersize=1.5, linestyle='-', linewidth=1, zorder=3)
        p_VM = ax.plot(t_sec, V_M, color='green', marker='o', markersize=1.5, linestyle='-', linewidth=1, zorder=2)
        p_VJM = ax.plot(t_sec, V_M - V_J, marker='o', markersize=1.5, linestyle='-', linewidth=1, zorder=2)
        p_VFM = ax.plot(t_sec, V_M - V_LG, marker='o', markersize=1.5, linestyle='-', linewidth=1, zorder=2)

        ax.legend(handles=[p_VJ[0], p_VM[0], p_VJM[0], p_VFM[0]], labels=['VJ', 'VM', 'VJM', 'VFM'], loc='upper left', bbox_to_anchor=(0.01, 1), frameon=True)
        ax.set_xlabel('t(s)', fontsize=10)
        ax.set_ylabel('V(mV)', fontsize=10)

        plt.show()