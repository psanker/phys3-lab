import numpy as np
from astropy import units as u
from astropy import constants as const
import sympy as sp

# relevant constants
PI     = np.pi
TWO_PI = 2*PI

R      = const.R
g      = const.g0

pa_si  = [u.Pa, u.N / (u.m)**2, lambda x: 1.0 * x, lambda x: 1.0 * x]
kpa_j  = [u.kPa * (u.m)**3, u.J, lambda x: 1e3 * x, lambda x: x / 1e3]

############################
# FUNCTIONS
############################
def vol(r, h):
    vol = (PI * (r[0] **2) * h[0])
    sv2 = (TWO_PI*r[0]*h[0] * r[1])**2 + (PI*l yea(r[0]**2) * h[1])**2

    return [vol, np.sqrt(sv2)]

def moles(r, h, T):
    V = vol(r, h)

    n   = (p_room * V[0]) / (R * T[0])
    sn2 = (((p_room/(R*T[0])) * V[1])**2) + (((-1*p_room*V[0])/(R*T[0]**2)) * T[1])**2
    return [n, np.sqrt(sn2)]

def work(points, h0):
    if len(points) is not 4:
        raise Exception('you dun goofed')

    a, b, c, d = points

    vol_a = vol(r_pis.to(u.m), [(h0[0] + a[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)])
    vol_b = vol(r_pis.to(u.m), [(h0[0] + b[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)])
    vol_c = vol(r_pis.to(u.m), [(h0[0] + c[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)])
    vol_d = vol(r_pis.to(u.m), [(h0[0] + d[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)])

    w1   = c[1] * (vol_c[0] - vol_b[0])
    sw12 = ((vol_c[1]**2 + vol_b[1]**2) * c[1]**2) + ((vol_c[0] - vol_b[0]) * (0.1 * u.kPa))**2

    w2   = a[1] * (vol_a[0] - vol_d[0])
    sw22 = ((vol_a[1]**2 + vol_d[1]**2) * a[1]**2) + ((vol_a[0] - vol_d[0]) * (0.1 * u.kPa))**2

    W    = w1 + w2
    sW2  = sw12 + sw22

    return [W.to(u.J, equivalencies=kpa_j), np.sqrt(sW2).to(u.J, equivalencies=kpa_j)]

def diff_vol(points, h0):
    if len(points) is not 4:
        raise Exception('you really can\'t count, m8')

    a, b, c, d = points

    vol_a = vol(r_pis.to(u.m), [(h0[0] + a[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)])
    vol_b = vol(r_pis.to(u.m), [(h0[0] + b[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)])
    vol_c = vol(r_pis.to(u.m), [(h0[0] + c[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)])
    vol_d = vol(r_pis.to(u.m), [(h0[0] + d[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)])

    diff0   = vol_c[0] - vol_b[0]
    s2diff0 = vol_c[1]**2 + vol_b[1]**2

    diff1   = vol_a[0] - vol_d[0]
    s2diff1 = vol_a[1]**2 + vol_d[1]**2

    return [diff0 + diff1, np.sqrt(s2diff0 + s2diff1)]


def gravity(points, h0, M):
    if len(points) is not 4:
        raise Exception('jeez dude')

    a, b, c, d = points

    h_a = [(h0[0] + a[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)]
    h_b = [(h0[0] + b[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)]
    h_c = [(h0[0] + c[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)]
    h_d = [(h0[0] + d[0]).to(u.m), (h0[1] + 1*u.mm).to(u.m)]

    del_h1   = h_b[0] - h_a[0]
    s2del_h1 = h_a[1]**2 + h_b[1]**2

    del_h2   = h_d[0] - h_c[0]
    s2del_h2 = h_c[1]**2 + h_d[1]**2

    Wi   = ((M + m_pis[0]).to(u.kg)) * g * del_h1
    sWi2 = (g * del_h1 * (m_pis[1]).to(u.kg))**2 + s2del_h1*((M + m_pis[0]).to(u.kg) * g)**2

    Wo   = (M + m_pis[0]).to(u.kg) * g * del_h2
    sWo2 = (g * del_h2 * (m_pis[1]).to(u.kg))**2 + s2del_h2*((M + m_pis[0]).to(u.kg) * g)**2

    return [(Wi + Wo).to(u.J), (np.sqrt(sWi2 + sWo2)).to(u.J)]

############################
# DATA
############################

r_pis  = np.array([16.25, 0.1]) * u.mm
m_pis  = np.array([35.0, 0.6]) * u.g

r_cham = np.array([2.000, 0.002]) * u.cm
h_cham = np.array([8.810, 0.002]) * u.cm

T_room = np.array([25.47 + 273.15, 0.01]) * u.K
p_room = 1.013e5 * (u.N / (u.m)**2)

# Coords in (x, p)
# For safety's sake, give all data points a σx of ±1mm considering the wobbliness of the terminal data points
# and give all points a σp of ± 0.1 kPa

# Trial 1
t1_r1_h0 = [35 * u.mm, 1*u.mm]
t1_r1_a  = [0 * u.m, .25 * u.kPa]
t1_r1_b  = [-0.003 * u.m, 2.7 * u.kPa]
t1_r1_c  = [0.012 * u.m, 2.7 * u.kPa]
t1_r1_d  = [0.015 * u.m, .25 * u.kPa]

t1_r2_h0 = [20 * u.mm, 1*u.mm]
t1_r2_a  = [0 * u.m, .25 * u.kPa]
t1_r2_b  = [-0.003 * u.m, 2.7 * u.kPa]
t1_r2_c  = [0.0105 * u.m, 2.7 * u.kPa]
t1_r2_d  = [0.013 * u.m, .25 * u.kPa]

t1_r3_h0 = [36 * u.mm, 1*u.mm]
t1_r3_a  = [0 * u.m, .25 * u.kPa]
t1_r3_b  = [-0.003 * u.m, 2.7 * u.kPa]
t1_r3_c  = [0.005 * u.m, 2.7 * u.kPa]
t1_r3_d  = [0.009 * u.m, .25 * u.kPa]

# Trial 2
t2_r1_h0 = [34 * u.mm, 1*u.mm]
t2_r1_a  = [0 * u.m, .20 * u.kPa]
t2_r1_b  = [-0.002 * u.m, 1.45 * u.kPa]
t2_r1_c  = [0.008 * u.m, 1.45 * u.kPa]
t2_r1_d  = [0.010 * u.m, .20 * u.kPa]

t2_r2_h0 = [35 * u.mm, 1*u.mm]
t2_r2_a  = [0 * u.m, .20 * u.kPa]
t2_r2_b  = [-0.002 * u.m, 1.5 * u.kPa]
t2_r2_c  = [0.006 * u.m, 1.5 * u.kPa]
t2_r2_d  = [0.007 * u.m, .20 * u.kPa]

t2_r3_h0 = [35 * u.mm, 1*u.mm]
t2_r3_a  = [0 * u.m, .20 * u.kPa]
t2_r3_b  = [-0.002 * u.m, 1.5 * u.kPa]
t2_r3_c  = [0.009 * u.m, 1.5 * u.kPa]
t2_r3_d  = [0.0105 * u.m, .20 * u.kPa]
