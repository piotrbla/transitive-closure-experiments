from islpy import *
from islplot.plotter import *
import matplotlib as mpl
import matplotlib.pyplot as plt

fig_poly = plt.figure()
domain = Set("{[i,j]: 0 <= i < 7 and i+1 <= j < 7 }")
#and 0 <= k <= -2 - i + j
# plot_set_points(Set("{S_0[i,j,k] : i >=0 and j < 7 and 0 <= k <= -2 - i + j}"),
#                marker=".", color="red", size=3)
# k = -2 - i + j
# plot_set_points(Set("{S_1[i, j] : i >= 0 and 2 + i <= j <= -2 + 7}"),
#                 marker="o", color="blue")
#dependences = Map("{[i,j]-> [i+1,j+1]; [i,j]-> [i+1,j]}")
#tiling = Map("{[i,j] -> [floor(i/2), floor(j/2)]}")
#space = Map("{[i,j] -> [i,i+j]}")
plot_domain(domain)
#plot_domain(domain, dependences, tiling, space)
#fig_poly.savefig("mcc.pdf", dpi=300)
fig_poly.savefig("mcc.png", dpi=300)