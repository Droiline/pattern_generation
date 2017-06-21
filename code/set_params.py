import reaction_equations as eq
from collections import OrderedDict

#this was getting too long and messy in the run_ah file
def set_params(test_key):
    model_type = 'basic'
    if test_key is 'meinhardt2_1':
        reaction_f = (eq.meinhardt2_1_a, eq.meinhardt2_1_i)
        # values to be explored for each parameter
        param_ranges = OrderedDict([
            ('dt', [1]),
            ('dx', [1]),
            ('innita',[6]),
            ('inniti',[4]),
            ('diffa',[0.3]),
            ('diffi',[0.02, 0.3, 0.4]),
            ('proda',[20]),
            ('prodi',[10]),
            ('rema',[5]),
            ('remi',[4]),
            ('sdens_0',[3]),
            ('sdens_var',[0]), #source density is uniform
            ('sdens',[0])
        ])
    elif test_key is 'meinhardt_paper1':
        reaction_f = (eq.meinhardt_paper1_a,eq.meinhardt_paper1_i)
        # values to be explored for each parameter
        param_ranges = OrderedDict([
            ('dt',[0]),
            ('innita',[1]),
            ('inniti',[1]),
            ('diffa',[0.002]),
            ('diffi',[0.4]),
            ('proda',[0.001]),
            ('prodi',[0.015]),
            ('rema',[0.01]),
            ('remi',[0]),
            ('sdens_0',[0.01]), # source density
            ('sdens_var',[0.025]),
            ('sdens',[0]), # empty, will be populated in model function.
            ('sata',[0.08]) # activator saturation
        ])
    elif test_key is 'meinhardt_paper2':
        reaction_f = (eq.meinhardt_paper2_a,eq.meinhardt_paper2_i)
        # values to be explored for each parameter
        param_ranges = OrderedDict([
            ('dt',[0]),
            ('innita',[0.1]),
            ('inniti',[0.1]),
            ('diffa',[0.1]),
            ('diffi',[0]),
            ('proda',[0.02]),
            ('prodi',[0.0075]),
            ('rema',[0.05]),
            ('remi',[0.03]),
            ('sdens_0',[0.1]),
            ('sdens_var',[0]),
            ('sdens',[0]), # empty, will be populated in model function.
            ('sata',[0.0004])
        ])
    elif test_key is 'meinhardt_branch':
        reaction_f = (eq.meinhardt_branch_a, eq.meinhardt_branch_i, eq.meinhardt_branch_r)
        model_type = 'globali'

        param_ranges = OrderedDict([
            ('dt',[5]),
            ('dx',[1]),
            ('innita',[0.5]),
            ('inniti',[0.5]),
            ('diffa',[0.015]),
            ('diffi',[0]),
            ('proda',[0.05, 0.1, 0.5]),
            ('prodi',[0.03, 0.06, 0.3]),
            ('rema',[0.05]),
            ('remi',[0.0075]),
            ('sdens_0',[0.02]),#rho_0 or r
            ('sdens_var',[0.14]),
            ('sdens',[0]),#rho
            ('sata',[0.004]), #kappa
            ('R',[0.05])
        ])
    elif test_key is 'meinhardt5_2':
        reaction_f = (eq.meinhardt5_2_a, eq.meinhardt5_2_i, eq.meinhardt5_2_i2)
        model_type = 'duali'
        # values to be explored for each parameter
        param_ranges = OrderedDict([
            ('dt',[4]),
            ('dx',[3]),
            ('innita',[1]),
            ('inniti',[1]),
            ('inniti2',[1]),
            ('diffa',[0.005]),
            ('diffi',[0]),
            ('diffi2',[0.3]),
            ('proda',[0]),
            ('prodi',[0.01]),
            ('prodi2',[0.02]),
            ('rema',[0.03]),
            ('remi',[0]),
            ('remi2',[0]),
            ('sdens_0', [0.05]),
            ('sdens_var', [0]),
            ('sdens', [0]),
            ('sata',[0.05])
        ])
    elif test_key is 'meinhardt5_3':
        reaction_f = (eq.meinhardt5_3_a, eq.meinhardt5_3_i, eq.meinhardt5_3_i2)
        model_type = 'duali'
        # values to be explored for each parameter
        param_ranges = OrderedDict([
            ('dt', [9]),
            ('innita',[1]),
            ('inniti',[1]),
            ('inniti2',[1]),
            ('diffa',[0.01]),
            ('diffi',[0.006]),
            ('diffi2',[0.4]),
            ('proda',[0.05]),
            ('prodi',[0]),
            ('prodi2',[0]),
            ('rema',[0.02]),
            ('remi',[0.002]),
            ('remi2',[0.01]),
            ('sdens_0', [0.008]),
            ('sdens_var', [0.05]),
            ('sdens', [0])
        ])
    elif test_key is 'meinhardt5_4':
        reaction_f = (eq.meinhardt5_4_a, eq.meinhardt5_4_i, eq.meinhardt5_4_i2)
        model_type = 'duali'
        # values to be explored for each parameter
        param_ranges = OrderedDict([
            ('dt',[6]),
            ('innita',[1]),
            ('inniti',[1]),
            ('inniti2',[1]),
            ('diffa',[0.01]),
            ('diffi',[0.005]),
            ('diffi2',[0.4]),
            ('proda',[0.08]),
            ('prodi',[0]),
            ('prodi2',[0]),
            ('rema',[0.01]),
            ('remi',[0.0015]),
            ('remi2',[0.015]),
            ('sdens_0',[0.05]),
            ('sdens_var',[0]),
            ('sdens',[0]), # empty, will be populated in model function.
            ('sati', [0.1]),
            ('sati2',[1])
        ])
    else:
        raise ValueError('Bad test key')

    return param_ranges, reaction_f, model_type
