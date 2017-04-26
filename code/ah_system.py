import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib
import sys
import os

def save_result(results, params, reaction_f):
    plt.close('all')
    titles = ['Activator levels','Inhibitor levels','Second Inhibitor levels']
    fig, axes = plt.subplots(1, len(results))

    for i, axis in enumerate(axes):
        axis.set_title(titles[i])
        axis.get_xaxis().set_visible(False)
        axis.get_yaxis().set_visible(False)
        axis.imshow(results[i], interpolation='none', cmap='YlOrBr')

    # save the result in a folder named after the function
    funcdir = reaction_f[0].__name__[:-2]
    # the filename is the parameters used
    filename = ''
    params.pop('sdens') # sdens is random, so not needed in filename.
    for k,v in params.items():
        filename = filename + k + str(v) + '_'

    filename = filename[:-1] + '.eps'
    dirpath = os.path.join('..','images',funcdir)

    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

    filepath = os.path.join(dirpath, filename)
    plt.savefig(filepath, bbox_inches='tight')

def ah_model(reaction_f, p, size, timesteps, debug=False):
    #create the model space
    activ = [[0]*size for _ in range(timesteps)]
    inhib = [[0]*size for _ in range(timesteps)]
    working_row = [[0,0]]*size
    sdens = p['sdens_0'] + np.random.rand(size)*p['sdens_var']

    stability = (p['dx']*p['dx']) / (2*p['dt'])

    if debug:
        print('dt: ', p['dt'])
        print('dx: ', p['dx'])
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

        save_result([activ, inhib], p, reaction_f)
        return True

#two inhibitor model
def ah2_model(reaction_f, p, size, timesteps, debug=False):
    #create the model space
    activ = [[0]*size for _ in range(timesteps)]
    inhib = [[0]*size for _ in range(timesteps)]
    inhib2 = [[0]*size for _ in range(timesteps)]
    working_row = [[0,0,0]]*size
    sdens = p['sdens_0'] + 2*p['sdens_var']*(np.random.rand(size)-0.5)

    stability = (p['dx']*p['dx']) / (2*p['dt'])

    if debug:
        print('dt: ', p['dt'])
        print('dx: ', p['dx'])
        print("diffusion rates should be <= ", stability)
        print("Diffa: ", p['diffa'])
        print("Diffi: ", p['diffi'])
        print("Diffi2: ", p['diffi2'])

    if p['diffa'] > stability or p['diffi'] > stability or p['diffi2'] > stability:
        raise ValueError("System is unstable, diffusion rates are " +
            str(p['diffa']) + ", " + str(p['diffi']) +
            " and " + str(p['diffi2']) +
            " and should be less than " + str(stability))
        return False
    else:
        #set initial values
        activ[0] = [p['innita']] * size
        inhib[0] = [p['inniti']] * size
        inhib2[0] = [p['inniti2']] * size

        #start iterating from after initial conditions
        for t in range(1,timesteps):
            for c in range(size):
                #deal with edge cases using the neutral flow method
                if c == 0:
                    old_activ_l = activ[t-1][c]
                    old_inhib_l = inhib[t-1][c]
                    old_inhib2_l = inhib2[t-1][c]
                else:
                    old_activ_l = activ[t-1][c-1]
                    old_inhib_l = inhib[t-1][c-1]
                    old_inhib2_l = inhib2[t-1][c-1]

                if c == size-1:
                    old_activ_r = activ[t-1][c]
                    old_inhib_r = inhib[t-1][c]
                    old_inhib2_r = inhib2[t-1][c]
                else:
                    old_activ_r = activ[t-1][c+1]
                    old_inhib_r = inhib[t-1][c+1]
                    old_inhib2_r = inhib2[t-1][c+1]

                #find the change in activator and inhibitor due to diffusion
                diffused_a = p['diffa'] * ((old_activ_l + old_activ_r - 2*activ[t-1][c]) / (p['dx']*p['dx']))
                diffused_i = p['diffi'] * ((old_inhib_l + old_inhib_r - 2*inhib[t-1][c]) / (p['dx']*p['dx']))
                diffused_i2 = p['diffi2'] * ((old_inhib2_l + old_inhib2_r - 2*inhib2[t-1][c]) / (p['dx']*p['dx']))

                #find the change in activator and inhibitor due to reaction
                p['sdens'] = sdens[c]
                reacted_a = reaction_f[0](activ[t-1][c], inhib[t-1][c], inhib2[t-1][c], p)
                reacted_i = reaction_f[1](activ[t-1][c], inhib[t-1][c], inhib2[t-1][c], p)
                reacted_i2 = reaction_f[2](activ[t-1][c], inhib[t-1][c], inhib2[t-1][c], p)
                #find new activator value
                activ[t][c] = activ[t-1][c] + p['dt'] * (diffused_a + reacted_a)
                inhib[t][c] = inhib[t-1][c] + p['dt'] * (diffused_i + reacted_i)
                inhib2[t][c] = inhib2[t-1][c] + p['dt'] * (diffused_i2 + reacted_i2)

        save_result([activ, inhib, inhib2], p, reaction_f)
        return True
