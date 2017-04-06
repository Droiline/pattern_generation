from ah_system import model_ah
from collections import OrderedDict
import time

debug = True

# the reaction equations
# rho is not used in the simpler equations, but they must all have the same
# parameter set for this to work.
def simpleMeinhardt_a(oldA, oldI, p, rho=0):
    dA = p['sdens'] * ((oldA * oldA) / oldI) - p['rema'] * oldA
    return dA

def simpleMeinhardt_i(oldA, oldI, p, rho=0):
    dI = p['sdens'] * oldA * oldA - p['remi'] * oldI + p['prodi']
    return dI

def meinhardt1_a(oldA, oldI, p, rho):
    aStar2 = ((oldA*oldA)/(1 + p['k'] * oldA*oldA)) + p['rema']
    dA = (rho * oldI * aStar2) - p['mu'] * oldA
    return dA

def meinhardt1_i(oldA, oldI, p, rho):
    aStar2 = ((oldA*oldA)/(1 + p['k'] * oldA*oldA)) + p['rema']
    dI = p['sdens'] - rho * oldI * aStar2 - p['nu'] * oldI
    return dI

def meinhardt2_a(oldA, oldI, p, rho):
    aStar2 = (oldA*oldA)/(1 + p['k'] * oldA*oldA)
    dA = (rho * (aStar2 + p['rho_0']))/oldI - p['mu'] * oldA
    return dA

def meinhardt2_i(oldA, oldI, p, rho):
    aStar2 = (oldA*oldA)/(1 + p['k'] * oldA*oldA);
    dI = rho * aStar2 - p['nu'] * oldI + p['remi'];
    return dI

size = 2000
time_steps = 2500

params = OrderedDict([
    ('innita',2),
    ('inniti',2),
    ('dt',0.5),
    ('dx',0.4),
    ('diffa',0.1),
    ('diffi',0.0),
    ('proda',6),
    ('prodi',4),
    ('rema',0.02),
    ('remi',0.0075),
    ('sdens',0.015),
    ('rho_0',0.02),
    ('rho_var',0.14),
    ('k',0.0004),
    ('mu',0.05),
    ('nu',0.03)
])

if debug:
    start = time.time()
try:
    model_ah(meinhardt2_a, meinhardt2_i, params, size, time_steps, debug)
except ValueError as e:
    print("ERROR: ", e)

if debug:
    end = time.time()
    print('time: ', end-start)
