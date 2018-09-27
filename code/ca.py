# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 10:19:38 2016

@author: droilin
"""

import sys
import numpy as np
from matplotlib import pyplot as plt

rule_number = 86
number_of_cells = 50
iterations = 50
random_init_conditions = False

rule_keys = ["111", "110", "101", "100",
         "011", "010", "001", "000"]

rule_string = "{0:b}".format(rule_number).zfill(8)

rules = dict(zip(rule_keys, rule_string))

if random_init_conditions:
    new_cells = [str(i) for i in list(np.random.choice(2, number_of_cells))]
else:
    new_cells = ['0'] * number_of_cells
    new_cells[len(new_cells)//2] = '1'

old_cells = list(new_cells)

output = np.zeros((iterations+1, number_of_cells))
output[0] = old_cells

for row in range(iterations):
    for i, cell in enumerate(new_cells):
        if (i == 0):
            rule = old_cells[len(old_cells)-1] + old_cells[i] + old_cells[i+1]
        elif (i == number_of_cells-1):
            rule = old_cells[i-1] + old_cells[i] + old_cells[0]
        else:
            rule = old_cells[i-1] + old_cells[i] + old_cells[i+1]

        cell = rules[rule]
        new_cells[i] = cell

    old_cells = list(new_cells)
    output[row+1] = old_cells

fig, ax = plt.subplots()
im = ax.pcolormesh(output, cmap='Greys')
ax.set_title('Cellular automata rule %i' % rule_number)
ax.axis('off')
plt.axis('equal')
plt.show()
