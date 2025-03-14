# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 23:52:14 2020

@author: Hong Zhu
"""

import numpy as np
import matplotlib.pyplot as plt


T1 = []
T2 = []
T3 = []
T4 = []
C = []
m = []
E = []
chi = []

with open("ising_model_C.txt", 'r') as fin_C:
    for words in fin_C:
        word = words.split()
        T1.append(float(word[0]))
        C.append(float(word[1]))

with open("ising_model_m.txt", 'r') as fin_m:
    for words in fin_m:
        word = words.split()
        T2.append(float(word[0]))
        m.append(float(word[1]))

with open("ising_model_E.txt", 'r') as fin_E:
    for words in fin_E:
        word = words.split()
        T3.append(float(word[0]))
        E.append(float(word[1]))

with open("ising_model_chi.txt", 'r') as fin_chi:
    for words in fin_chi:
        word = words.split()
        T4.append(float(word[0]))
        chi.append(float(word[1]))

plt.figure(figsize=(4,8))
plt.subplot(411)
plt.plot(T3, E, 'r.', label='C')
plt.xlim(0, 4)
plt.ylabel(r"$\langle E \rangle$")

plt.gcf().subplots_adjust(left=0.2)
plt.subplot(412)
plt.plot(T1, C, 'r.', label='C')
plt.xlim(0, 4)
plt.ylabel(r"$C$")

plt.subplot(413)
plt.plot(T2, m, 'r.', label='C')
plt.xlim(0, 4)
plt.ylabel(r"$\langle |m| \rangle$")

plt.subplot(414)
plt.plot(T4, chi, 'r.', label=r'$\chi$')
plt.xlim(0, 4)
plt.xlabel(r"$T$")
plt.ylabel(r"$\chi$")

plt.savefig("75x75.png", dpi=150, bbox_inches='tight')
plt.show()


