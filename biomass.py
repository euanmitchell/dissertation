#####################################
# Impact of ground elevation bias on
# biomass using MCH relationship.
#####################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Set plot parameters
plt.rcParams['figure.figsize']=(8,6)
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams['axes.labelsize']=18
plt.rcParams['axes.labelpad']=8.0
plt.rcParams['axes.titlesize']=20
#print(plt.rcParams)

# Calculate biomass (ACD/AGBD) as a function of mean canopy height for temperate and tropical forests
mch=np.arange(5,71,1)
acdTemp=0.378*(mch**2)
acdTrop=0.844*(mch**2)

# Calculate mch with bias
mchTempBias=mch+1
mchTropBias=mch-4
tropCorr=mch-1.43

# Calculate biomass with bias
acdTempBias=0.378*(mchTempBias**2)
acdTempBiasPct=(acdTempBias/acdTemp)*100
acdTropBias=0.844*(mchTropBias**2)
acdTropBiasPct=(acdTropBias/acdTrop)*100
tropCorrBias=0.844*(tropCorr**2)
tropCorrBiasPct=(tropCorrBias/acdTrop)*100

fig,ax=plt.subplots()
ax.add_patch(Rectangle((5,80),65,40,facecolor='lightgrey'))
ax.plot(mch,acdTempBiasPct,label='Temperate Forest 1 m Bias')
ax.plot(mch,acdTropBiasPct,label='Tropical Forest -4 m Bias')
ax.plot(mch,tropCorrBiasPct,label='Tropical Forest -1.43 m Bias',linestyle='--',color='orange')
ax.axhline(y=100,color='black',linestyle='--')
plt.legend()
plt.xlim([5,70])
plt.ylim([0,150])
plt.xlabel('Mean Canopy Height (m)')
plt.ylabel('Biomass (% of Reference Value)')
plt.savefig('../data/biomassBiasCorr.png',dpi=300)
plt.close()
plt.clf()
#plt.show()
