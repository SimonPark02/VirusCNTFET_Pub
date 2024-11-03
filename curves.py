"""Define special functions"""


from scipy.special import expit
from numpy import absolute
from numpy import heaviside
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