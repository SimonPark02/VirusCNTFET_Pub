"""
Dimesion                 Unit
--------------------     --------
area                     cm^2
volume                   cm^3
electric potential       mV
time                     ms
specific conductance     mS/cm^2
specific capacitance     uF/cm^2
ion concentration        mM
"""

from scipy.constants import pi

# Total membrane area
AM = (pi * 30 * 300 + 2 * pi * 225) * 10 ** (-14)

# Ratio between attached & total membrane area
r = 0.3

# Attached membrane area
AJ = r * AM

# Specific capacitance of free membrane
cm = 1.0

# Specific capacitance of juction area
cJ = 2.4

# Specific conductance of junction area
gJ = 0.1

# Specific conductance of attached membrane
gKJ = 28.0

# Specific conductance of free membrane
gKF = 8.4

# Channel activation parameters
Vh = -45.6
km = 3.0

# Resting membrane potential
Vrest = -36

# Potassiun reversal potential
VK = -48