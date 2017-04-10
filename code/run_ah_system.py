from ah_system import model_ah
from collections import OrderedDict
import time
import itertools as it
import numpy as np

debug = False
full_run = False
test_key = 'rmeinhardt2'

# the reaction equations
def meinhardt_chap2_a(oldA, oldI, p):
    dA = p['sdens'] * (((oldA * oldA) / oldI) + p['proda']) - p['rema'] * oldA
    return dA

def meinhardt_chap2_i(oldA, oldI, p):
    dI = p['sdens'] * oldA * oldA - p['remi'] * oldI + p['prodi']
    return dI

def rmeinhardt1_a(oldA, oldI, p):
    aStar2 = ((oldA*oldA)/(1 + p['k'] * oldA*oldA)) + p['rema']
    dA = (['rho'] * oldI * aStar2) - p['mu'] * oldA
    return dA

def rmeinhardt1_i(oldA, oldI, p):
    aStar2 = ((oldA*oldA)/(1 + p['k'] * oldA*oldA)) + p['rema']
    dI = p['sdens'] - ['rho'] * oldI * aStar2 - p['nu'] * oldI
    return dI

def rmeinhardt2_a(oldA, oldI, p):
    aStar2 = (oldA*oldA)/(1 + p['k'] * oldA*oldA)
    dA = (p['rho'] * (aStar2 + p['rho_0']))/oldI - p['mu'] * oldA
    return dA

def rmeinhardt2_i(oldA, oldI, p):
    aStar2 = (oldA*oldA)/(1 + p['k'] * oldA*oldA);
    dI = p['rho'] * aStar2 - p['nu'] * oldI + p['remi'];
    return dI

if test_key is 'meinhardt_chap2':
    # values to be explored for each parameter
    param_ranges = OrderedDict([
        ('dt',[1]),
        ('dx',[1]),
        ('innita',np.arange(1,10,2)),
        ('inniti',np.arange(1,10,2)),
        ('diffa',np.arange(0,1,0.1)),
        ('diffi',np.arange(0,1,0.1)),
        ('proda',np.arange(1, 10, 1)),
        ('prodi',np.arange(1, 10, 1)),
        ('rema',np.arange(0.01,0.1,0.01)),
        ('remi',np.arange(0.01,0.1,0.01)),
        ('sdens',np.arange(0.01,0.1,0.01))
    ])
    reaction_f = (meinhardt_chap2_a, meinhardt_chap2_i)
elif test_key is 'rmeinhardt1':
    # values to be explored for each parameter
    param_ranges = OrderedDict([
        ('dt',[1]),
        ('dx',[1]),
        ('innita',np.arange(1,10,2)),
        ('inniti',np.arange(1,10,2)),
        ('diffa',np.arange(0,1,0.1)),
        ('diffi',np.arange(0,1,0.1)),
        ('proda',np.arange(1, 10, 1)),
        ('prodi',np.arange(1, 10, 1)),
        ('rema',np.arange(0.01,0.1,0.01)),
        ('remi',np.arange(0.01,0.1,0.01)),
        ('sdens',np.arange(0.01,0.1,0.01)),
        ('rho_0',np.arange(0.01,0.1,0.01)),
        ('rho_var',np.arange(0.1,0.2,0.01)),
        ('rho',[0]), # empty, will be populated in model function.
        ('k',np.arange(0.0001,0.001,0.0001)),
        ('mu',np.arange(0.01,0.1,0.01)),
        ('nu',np.arange(0.01,0.1,0.01))
    ])
    reaction_f = (rmeinhardt1_a,rmeinhardt1_i)
else:
    # values to be explored for each parameter
    param_ranges = OrderedDict([
        ('dt',[1]),
        ('dx',[1]),
        ('innita',np.arange(1,10,2)),
        ('inniti',np.arange(1,10,2)),
        ('diffa',np.arange(0,1,0.1)),
        ('diffi',np.arange(0,1,0.1)),
        ('proda',np.arange(1, 10, 1)),
        ('prodi',np.arange(1, 10, 1)),
        ('rema',np.arange(0.01,0.1,0.01)),
        ('remi',np.arange(0.01,0.1,0.01)),
        ('sdens',np.arange(0.01,0.1,0.01)),
        ('rho_0',np.arange(0.01,0.1,0.01)),
        ('rho_var',np.arange(0.1,0.2,0.01)),
        ('rho',[0]), # empty, will be populated in model function.
        ('k',np.arange(0.0001,0.001,0.0001)),
        ('mu',np.arange(0.01,0.1,0.01)),
        ('nu',np.arange(0.01,0.1,0.01))
    ])
    reaction_f = (rmeinhardt2_a,rmeinhardt2_i)

# in theory this generator will take a list of lists of possible
# parameter values and yield every possible combination of params.
# in theory.
# param_list is a list of tuples and param_set is an OrderedDict.
# note that param_list does not need to be in the same order as
# param_set.
def product(param_list, param_set):
    for val in param_list[0][1]:
        # set the value of the first parameter
        param_set[param_list[0][0]] = val
        # if this isn't the last parameter, recurse
        if len(param_list) is not 1:
            for prod in product(param_list[1:], param_set):
                yield prod
        # else we have a complete parameter set and can yield it
        else:
            yield param_set

def set_steps(diffa, diffi):
    # an attempt at automatically setting the correct dx and dt
    # depending on the diffusion rates.
    # if we fix dx = 1 then dt < 1/2*diff
    dx = 1
    diff = max(diffa, diffi) if max(diffa, diffi) > 0 else 1
    stable = 0.9 * (1/(2*diff))
    dt = stable if stable < 1 else 1

    # alternativeley, if we fix dx = dt then dt > 2 * diff
    return dx, dt

size = 2000
time_steps = 2500

# list all parameter combinations
#param_sets = list(it.product(*param_ranges.values()))
# this didn't work: python replied 'Killed', wasn't beefy enough
# going to try using many loops instead.
# no, many loops is ugly. recursion?
# recursive generator works nicely.

empty_p = OrderedDict(zip(param_ranges.keys(),[0]*len(param_ranges.keys())))

start = time.time()

if full_run:
    for params in product(list(param_ranges.items()), empty_p):
        params['dx'], params['dt'] = set_steps(params['diffa'],params['diffi'])

        try:
            model_ah(reaction_f, params, size, time_steps, debug)
        except ValueError as e:
            print("ERROR: ", e)
else:
    params = [2,2,1,1,0.1,0,6,4,0.02,0.0075,0.015,0.02,0.14,0,0.0004,0.05,0.03]
    params = OrderedDict(zip(param_ranges.keys(),params))
    params['dx'], params['dt'] = set_steps(params['diffa'],params['diffi'])

    try:
        model_ah(reaction_f, params, size, time_steps, debug)
    except ValueError as e:
        print("ERROR: ", e)

end = time.time()
print('time: ', end-start)

#params = OrderedDict([
#    ('innita',2),
#    ('inniti',2),
#    ('dt',0.5),
#    ('dx',0.4),
#    ('diffa',0.1),
#    ('diffi',0.0),
#    ('proda',6),
#    ('prodi',4),
#    ('rema',0.02),
#    ('remi',0.0075),
#    ('sdens',0.015),
#    ('rho_0',0.02),
#    ('rho_var',0.14),
#    ('k',0.0004),
#    ('mu',0.05),
#    ('nu',0.03)
#])
