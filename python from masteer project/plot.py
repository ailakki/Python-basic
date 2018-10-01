# if seaborn install works
import seaborn as sns
import matplotlib.pyplot as plt
import pylab
import sys

sns.set_style("whitegrid")
fil = raw_input("File name of genome to plot :")
f = open(fil,"r")

name = "None"

fig, ax =plt.subplots(3,2)
#axes = axes.flatten()
n=0
name= 'a'
for line in f:
	dat=[]
	if line.startswith('#'): 
		name = line[1:]
	else:
		#print line
		l = line.strip() 
		d = l.split()
		for x in d:
			dat.append(float(x))
		ax_curr = ax[abs(n/2),n%2]
		ax_curr.set_title(name)
		sns.violinplot(data=[dat],orient="h", color='b', ax=ax_curr)
		plt.xlim(-100,100)
		n = n+1
#fig.subplots_adjust(hspace=0.3) 
plt.show()
