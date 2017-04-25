import reaction_equations as eq
from collections import OrderedDict

#this was getting too long and messy in the run_ah file
def set_params(test_key, full_run=True):
    model_type = 'basic'
    if test_key is 'meinhardt2_1':
        reaction_f = (eq.meinhardt2_1_a, eq.meinhardt2_1_i)
        # values to be explored for each parameter
        param_ranges = OrderedDict([
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
    elif test_key is 'meinhardt_paper1':
        reaction_f = (eq.meinhardt_paper1_a,eq.meinhardt_paper1_i)
        # values to be explored for each parameter
        if full_run:
            param_ranges = OrderedDict([
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
        else:
            param_ranges = OrderedDict([
                ('innita',[2]),
                ('inniti',[2]),
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
        if full_run:
            param_ranges = OrderedDict([
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
                ('sata',[0.0004])
            ])
        else:
            param_ranges = OrderedDict([
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
                ('sata',[0.0004])
            ])
    elif test_key is 'meinhardt5_3':
        reaction_f = (eq.meinhardt5_3_a, eq.meinhardt5_3_i, eq.meinhardt5_3_i2)
        model_type = 'duali'
        # values to be explored for each parameter
        if full_run:
            param_ranges = OrderedDict([
                ('dx', [0.3]),
                ('dt', [3]),
                ('innita',[0.01,0.1,1,10]),
                ('inniti',[0.01,0.1,1]),
                ('inniti2',[0.01,0.1,1]),
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
                ('sdens_var', [0,0.00015]),
                ('sdens', [0])
            ])
        else:
            param_ranges = OrderedDict([
                ('innita',[0.1]),
                ('inniti',[0.1]),
                ('inniti2',[0.1]),
                ('diffa',[0.01]),
                ('diffi',[0.006]),
                ('diffi2',[0.4]),
                ('proda',[0.05]),
                ('prodi',[0]),
                ('prodi2',[0]),
                ('rema',[0.02]),
                ('remi',[0.002]),
                ('remi2',[0.01])
            ])
    elif test_key is 'meinhardt5_4':
        reaction_f = (eq.meinhardt5_4_a, eq.meinhardt5_4_i, eq.meinhardt5_4_i2)
        model_type = 'duali'
        # values to be explored for each parameter
        if full_run:
            param_ranges = OrderedDict([
                ('dx',[3])
                ('innita',[0.1]),
                ('inniti',[0.1]),
                ('inniti2',[0.1]),
                ('diffa',[0.01]),
                ('diffi',[0.006]),
                ('diffi2',[0.4]),
                ('proda',[0.05]),
                ('prodi',[0]),
                ('prodi2',[0]),
                ('rema',[0.02]),
                ('remi',[0.002]),
                ('remi2',[0.01]),
                ('sdens_0',[0.05]),
                ('sdens_var',[0.015]),
                ('sdens',[0]), # empty, will be populated in model function.
                ('sata',[0.0004]),
                ('sati', [0.0004]),
                ('sati2',[0.0004])
            ])
        else:
            param_ranges = OrderedDict([
                ('innita',[0.1]),
                ('inniti',[0.1]),
                ('inniti2',[0.1]),
                ('diffa',[0.1]),
                ('diffi',[0]),
                ('diffi2',[0]),
                ('proda',[0.02]),
                ('prodi',[0.0075]),
                ('prodi2',[0.0075]),
                ('rema',[0.05]),
                ('remi',[0.03]),
                ('remi2',[0.03]),
                ('sdens_0',[0.05]),
                ('sdens_var',[0.015]),
                ('sdens',[0]), # empty, will be populated in model function.
                ('sata',[0.0004])
            ])
    else:
        raise ValueError('Bad test key')

    return param_ranges, reaction_f, model_type
