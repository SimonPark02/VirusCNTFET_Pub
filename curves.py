"""Define special functions"""


from scipy.special import expit
from numpy import absolute
from numpy import heaviside
from typing import Callable
from parameters import *


def m_inf(V, Vh, k):
    return expit((V - Vh) / k)


def dep(w, x):
    return 1 / (2 * w) * (absolute(x) - absolute(x - w) + w)


def ddep(w, x):
    return 1 / w * (heaviside(x, 0.5) - heaviside(x - w, 0.5))


def tri(w, x):
    return 1 / w * (0.5 * (absolute(x - w) + absolute(x + w)) - absolute(x))


def dtri(w, x):
    return 1 / w * (heaviside(x - w, 0.5) + heaviside(x + w, 0.5) - 2 * heaviside(x, 0.5))


#TODO: Add documentation
def create_step_voltage(
    v_step: float,
    v_sd: float,
    t_polarization: float,
    polarization_delay: float
) -> tuple[Callable, Callable, Callable]:
    """Creates voltage functions for step voltage"""
    def Vlg(t):
        """Liquid-gate voltage"""
        return v_step * dep(polarization_delay, t - (t_polarization - polarization_delay))

    def dVlgdt(t):
        """Time derivative of liquid-gate voltage"""
        return v_step * ddep(polarization_delay, t - (t_polarization - polarization_delay))

    def dVsdt(t):
        """Time derivative of source-drain voltage"""
        return v_sd * ddep(polarization_delay, t - (t_polarization - polarization_delay))

    return Vlg, dVlgdt, dVsdt


#TODO: Add documentation
def create_voltage_sweep(
    v_low: float,
    v_high: float,
    slope: float,
    v_sd: float,
    t_polarization: float,
    polarization_delay: float
) -> tuple[Callable, Callable, Callable]:
    """Creates voltage functions for cyclic voltametry"""
    _height = v_high - v_low
    _width = _height / slope

    def Vlg(t):
        """Liquid-gate voltage"""
        return v_low * dep(polarization_delay, t - (t_polarization - polarization_delay)) + _height * tri(_width, t - _width - t_polarization)

    def dVlgdt(t):
        """Time derivative of liquid-gate voltage"""
        return v_low * ddep(polarization_delay, t - (t_polarization - polarization_delay)) + _height * dtri(_width, t - _width - t_polarization)

    def dVsdt(t):
        """Time derivative of source-drain voltage"""
        return v_sd * ddep(polarization_delay, t - (t_polarization - polarization_delay))

    return Vlg, dVlgdt, dVsdt