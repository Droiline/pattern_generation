# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 10:19:38 2016

@author: droilin
"""

import sys
import numpy as np

if len(sys.argv) == 0:
    print("For help run: python ca.py -h")

rule_number = 86
number_of_cells = 25
iterations = 30
random_init_conditions = True

rule_keys = ["111", "110", "101", "100", 
         "011", "010", "001", "000"]

rule_number = "{0:b}".format(rule_number)
rule_number = str(rule_number).zfill(8)

rules = dict(zip(rule_keys, rule_number))

print(rule_number)
print(rules)

if random_init_conditions:
    new_cells = [str(i) for i in list(np.random.choice(2, number_of_cells))]
else:
    new_cells = [0] * number_of_cells
    new_cells[len(new_cells)//2] = 1

old_cells = list(new_cells)
print("s :", old_cells)

for row in range(iterations):
    print(row, ":", end=' ')
    for i, cell in enumerate(new_cells):
        if (i == 0):
            rule = old_cells[len(old_cells)-1] + old_cells[i] + old_cells[i+1]
        elif (i == number_of_cells-1):
            rule = old_cells[i-1] + old_cells[i] + old_cells[0]
        else:
            rule = old_cells[i-1] + old_cells[i] + old_cells[i+1]
        
        cell = rules[rule]
        new_cells[i] = cell
        
    print(new_cells)
    old_cells = list(new_cells)