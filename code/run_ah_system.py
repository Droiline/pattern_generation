from ah_system import model_ah
from collections import OrderedDict
import time
import itertools as it
import numpy as np

debug = False
full_run = False
test_key = 'meinhardt2'
size = 500
time_steps = 700

# the reaction equations
def meinhardt2_1_a(oldA, oldI, p):
    dA = p['sdens_0'] * (((oldA * oldA) / oldI) + p['proda']) - p['rema'] * oldA
    return dA

def meinhardt2_1_i(oldA, oldI, p):
    dI = p['sdens_0'] * oldA * oldA - p['remi'] * oldI + p['prodi']
    return dI

#an activator-substrate systems
def meinhardt_paper1_a(oldA, oldI, p):
    aStar2 = ((oldA*oldA)/(1 + p['sat'] * oldA*oldA)) + p['proda']
    dA = (p['sdens'] * oldI * aStar2) - p['rema'] * oldA
    return dA

def meinhardt_paper1_i(oldA, oldI, p):
    aStar2 = ((oldA*oldA)/(1 + p['sat'] * oldA*oldA)) + p['proda']
    dI = p['prodi'] - p['sdens'] * oldI * aStar2 - p['remi'] * oldI
    return dI

#an activator-inhibitor system
def meinhardt_paper2_a(oldA, oldI, p):
    aStar2 = (oldA*oldA)/(1 + p['sat'] * oldA*oldA)
    dA = (p['sdens'] * (aStar2 + p['proda']))/oldI - p['rema'] * oldA
    return dA

def meinhardt_paper2_i(oldA, oldI, p):
    aStar2 = (oldA*oldA)/(1 + p['sat'] * oldA*oldA);
    dI = p['sdens'] * aStar2 - p['remi'] * oldI + p['prodi'];
    return dI

if test_key is 'meinhardt2_1':
    # values to be explored for each parameter
    param_ranges = OrderedDict([
        ('dt',[1]),
        ('dx',[1]),
        ('innita',np.arange(1,10,2)),
        ('inniti',np.arange(1,10,2)),
        ('diffa',np.arange(0.01,2,0.2)),
        ('diffi',np.arange(0.0,2,0.2)),
        ('proda',np.arange(0.0, 100, 10)),
        ('prodi',np.arange(0, 100, 10)),
        ('rema',np.arange(0,10,2)),
        ('remi',np.arange(0,10,2)),
        ('sdens_0',np.arange(0,10,2)),
        ('sdens_var',[0]), #source density is uniform
        ('sdens',[0])
    ])
    reaction_f = (meinhardt2_1_a, meinhardt2_1_i)
elif test_key is 'meinhardt_paper1':
    reaction_f = (meinhardt_paper1_a,meinhardt_paper1_i)
    # values to be explored for each parameter
    if full_run:
        param_ranges = OrderedDict([
            ('dt',[1]),
            ('dx',[1]),
            ('innita',[1]),
            ('inniti',[1]),
            ('diffa',[0.002]),
            ('diffi',[0.4]),
            ('proda',[0.001]),
            ('prodi',[0.015]),
            ('rema',[0.01]),
            ('remi',[0]),
            ('sdens_0',[0.01]), # source density?
            ('sdens_var',[0.025]),
            ('sdens',[0]), # empty, will be populated in model function.
            ('sat',[0.08]) # activator saturation
            ])
    else:
        param_ranges = OrderedDict([
            ('dt',[1]),
            ('dx',[1]),
            ('innita',[2]),
            ('inniti',[2]),
            ('diffa',[0.002]),
            ('diffi',[0.4]),
            ('proda',[0.001]),
            ('prodi',[0.015]),
            ('rema',[0.01]),
            ('remi',[0]),
            ('sdens_0',[0.01]), # source density?
            ('sdens_var',[0.025]),
            ('sdens',[0]), # empty, will be populated in model function.
            ('sat',[0.08]) # activator saturation
        ])
else:
    reaction_f = (meinhardt_paper2_a,meinhardt_paper2_i)
    # values to be explored for each parameter
    if full_run:
        param_ranges = OrderedDict([
            ('dt',[1]),
            ('dx',[1]),
            ('innita',[0.1]),
            ('inniti',[0.1]),
            ('diffa',[0.1]),
            ('diffi',[0]),
            ('proda',[0.02]),
            ('prodi',[0.0075]),
            ('rema',[0.05]),
            ('remi',[0.03]),
            ('sdens_0',[0.05]),
            ('sdens_var',[0.015]),
            ('sdens',[0]), # empty, will be populated in model function.
            ('sat',[0.0004])
        ])
    else:
        param_ranges = OrderedDict([
            ('dt',[1]),
            ('dx',[1]),
            ('innita',[2]),
            ('inniti',[2]),
            ('diffa',[0.1]),
            ('diffi',[0]),
            ('proda',[0.02]),
            ('prodi',[0.0075]),
            ('rema',[0.05]),
            ('remi',[0.03]),
            ('sdens_0',[0.05]),
            ('sdens_var',[0.015]),
            ('sdens',[0]), # empty, will be populated in model function.
            ('sat',[0.0004])
        ])

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
#    dx = 1
#    diff = max(diffa, diffi) if max(diffa, diffi) > 0 else 1
#    stable = 0.9 * (1/(2*diff))
#    dt = stable if stable < 10 else 10

    # alternatively, if we fix dx = dt then dt > 2 * diff
    diff = min(diffa, diffi) if min(diffa, diffi) > 0 else 1
    stable = 1.1 * 2 * diff
    dx = stable
    dt = stable
    return dx, dt

runs = 0

# list all parameter combinations
#param_sets = list(it.product(*param_ranges.values()))
# this didn't work: python replied 'Killed', wasn't beefy enough
# going to try using many loops instead.
# no, many loops is ugly. recursion?
# recursive generator works nicely.

empty_p = OrderedDict(zip(param_ranges.keys(),[0]*len(param_ranges.keys())))

start = time.time()

for params in product(list(param_ranges.items()), empty_p):
    params['dx'], params['dt'] = set_steps(params['diffa'],params['diffi'])

    try:
        runs = runs + 1
        model_ah(reaction_f, params, size, time_steps, debug)
    except ValueError as e:
        runs = runs - 1
        print("ERROR: ", e)

end = time.time()
print('number of runs: ', runs)
print('time taken: ', end-start)
print('time per run: ', (end-start)/runs)

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
#    ('sat',0.0004),
#    ('mu',0.05),
#    ('nu',0.03)
#])
