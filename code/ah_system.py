import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib
import sys
import os

def model_ah(reaction_f, p, size, timesteps, debug=False):
    #create the model space
    activ = [[0]*size for _ in range(timesteps)]
    inhib = [[0]*size for _ in range(timesteps)]
    working_row = [[0,0]]*size
    sdens = p['sdens_0'] + np.random.rand(size)*p['sdens_var']

    stability = (p['dx']*p['dx']) / (2*p['dt'])

    if debug:
        print("diffusion rates should be <= ", stability)
        print("Diffa: ", p['diffa'])
        print("Diffi: ", p['diffi'])

    if p['diffa'] >= stability or p['diffi'] >= stability:
        raise ValueError("System is unstable, diffusion rates should be less than " + str(stability))
        return False
    else:
        #set initial values
        activ[0] = [p['innita']] * size
        inhib[0] = [p['inniti']] * size
#        activ[0][int(size/2)] = p['innita'] + 1
#        activ[0][int(size/2)+1] = p['innita'] + 1
#        inhib[0][int(size/2)] = p['inniti'] + 1
#        inhib[0][int(size/2)+1] = p['inniti'] + 1

        #start iterating from after initial conditions
        for t in range(1,timesteps):
            for c in range(size):
                #deal with edge cases using the neutral flow method
                if c == 0:
                    old_activ_l = activ[t-1][c]
                    old_inhib_l = inhib[t-1][c]
                else:
                    old_activ_l = activ[t-1][c-1]
                    old_inhib_l = inhib[t-1][c-1]

                if c == size-1:
                    old_activ_r = activ[t-1][c]
                    old_inhib_r = inhib[t-1][c]
                else:
                    old_activ_r = activ[t-1][c+1]
                    old_inhib_r = inhib[t-1][c+1]

                #find the change in activator and inhibitor due to diffusion
                diffused_a = p['diffa'] * ((old_activ_l + old_activ_r - 2*activ[t-1][c]) / (p['dx']*p['dx']))
                diffused_i = p['diffi'] * ((old_inhib_l + old_inhib_r - 2*inhib[t-1][c]) / (p['dx']*p['dx']))

                #find the change in activator and inhibitor due to reaction
                p['sdens'] = sdens[c]
                reacted_a = reaction_f[0](activ[t-1][c], inhib[t-1][c], p)
                reacted_i = reaction_f[1](activ[t-1][c], inhib[t-1][c], p)
                #find new activator value
                activ[t][c] = activ[t-1][c] + p['dt'] * (diffused_a + reacted_a)
                inhib[t][c] = inhib[t-1][c] + p['dt'] * (diffused_i + reacted_i)

        plt.close('all')
        fig, axes = plt.subplots(1, 2)
        axes[0].set_title('Activator levels')
        axes[0].get_xaxis().set_visible(False)
        axes[0].get_yaxis().set_visible(False)
        im1 = axes[0].imshow(activ, interpolation='none', cmap='YlOrBr')
        axes[1].set_title('Inhibitor levels')
        axes[1].get_xaxis().set_visible(False)
        axes[1].get_yaxis().set_visible(False)
        im2 = axes[1].imshow(inhib, interpolation='none', cmap='YlOrBr')

        # save the result in a folder named after the function
        funcdir = reaction_f[0].__name__[:-2]
        # the filename is the parameters used
        filename = ''
        p.pop('sdens') # sdens is random, so not needed.
        for k,v in p.items():
            filename = filename + k + str(v) + '_'

        filename = filename[:-1] + '.eps'
        dirpath = os.path.join('..','images',funcdir)

        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

        filepath = os.path.join(dirpath, filename)
        plt.savefig(filepath, bbox_inches='tight')

        return True
