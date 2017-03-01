# Whatever the heading is I put for these things

debug = False
# Given parameters
g      = 3*9.81  # [-]
L      = 70.0    # [m]
Lf1    = 5.0     # [m]
Lf2    = 31.2    # [m]
Lf3    = 11.0    # [m]
R      = 3.1     # [m]
hf     = 2.0     # [m]
ts     = 0.004   # [m]
tf     = 0.025   # [m]
tst    = 0.0012  # [m]
hst    = 0.015   # [m]
wst    = 0.02    # [m]
ns     = 36      # [-]
dtailz = 7.5     # [m]
dlgy   = 3       # [m]
Sx     = 6.1e5   # [m]
M      = 250000  # [kg]

# Calculated values
W = M*g  # weight
q = W / L  # q is assumed to be a 1-dimensional line distributed load
q_alt = W / (L*2*R)  # q_alt is a 2-dimensional area distributed load
Sy2 = W*(L/2 - Lf1) / Lf2
Sy1 = W - Sy2
M3 = Sx * dtailz
Sx2 = -1376410.256
Sx1 = 766410.2564
Mlg1 = 2299230.769
Mlg2 = -4129230.769 + 6405000
Mtail = -4575000

Lf02 = L - Lf1 - Lf2  # distance between z=0 and rear landing gear
Lf01 = L - Lf1  # distance between z=0 and front landing gear

# alt_version
external_forces_y = [-q, Sy1, Sy2]
external_locations_y = [(0., L), Lf01, Lf02] # this is the z-location for the forces in y

assert len(external_forces_y) == len(external_locations_y), 'Load values/locations are not of equal length.'
