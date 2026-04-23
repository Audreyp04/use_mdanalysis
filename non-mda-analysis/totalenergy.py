#Written by Audrey D. Prendergast
#Fully functional as of 4/22/2026

import matplotlib.pyplot as plt
import pyedr

#First must have a concatenated EDR file 
#gmx eneconv -f *.edr -o combined.edr 

path = 'combined.edr'
dic = pyedr.edr_to_dict(path)

print(dic)

time=dic['Time']/1000
potential=dic['Total Energy']

plt.plot(time,potential, c='orchid',ls='-')
plt.xlim(0,2000)
plt.xlabel('Time (ns)')
plt.ylabel('Total Energy (kJ/mol)')
plt.savefig('energy.png')
plt.close()