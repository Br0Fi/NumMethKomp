# Importanweisungen
import matplotlib.patches as mpatches
import numpy as np
import scipy as sci
import scipy.fftpack
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.axes as axes
from matplotlib import colors as mcolors
import math

# Konstanten fuer einheitliche Darstellung
fig_size = (10, 6)
fig_legendsize = 15
fig_labelsize = 15
matplotlib.rcParams.update({'font.size': fig_labelsize})
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

#Aufgabenspezifische Konstanten:
a = 2
b = 0.02
c = 0.0002
d = 0.8

#TODO:
# Ergebnis interpretieren -> vlt. doch im Phasenraum plotten -> siehe VL-Kram

# Hauptprogramm:


fig = plt.figure(figsize=fig_size)

plt.grid()
plt.legend(prop={'size':fig_legendsize})
plt.ylabel('x')
plt.xlabel("t")
plt.tick_params(labelsize=fig_labelsize)
# plt.savefig("Zettel1/figures/A1-bc.png")
plt.show()
# %%
