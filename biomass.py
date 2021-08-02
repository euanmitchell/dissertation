#####################################
# Impact of ground elevation bias on
# biomass using MCH relationship.
#####################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse

#######################################

# Defining the command line reading function
def readCommands():
    '''
    Read the arguments passed from the command line
    '''
    p=argparse.ArgumentParser(description=('Specify input ALS and GEDI data files and programme control'),
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--tempBias', dest='tempBias', type=float, default='', help=('The bias for temperate forests'))
    p.add_argument('--tempCorr', dest='tempCorr', type=float, default=0.0, help=('The corrected bias for temperate forests'))
    p.add_argument('--tropBias', dest='tropBias', type=float, default='', help=('The bias for tropical forests'))
    p.add_argument('--tropCorr', dest='tropCorr', type=float, default=0.0, help=('The corrected bias for tropical forests'))
    p.add_argument('--output', dest='output', type=str, default='',
        help=('Output path and filename for plot.\nDefault is not set'))
    cmdargs=p.parse_args()
    return cmdargs

args=readCommands()

#######################################

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
mchTempBias=mch+args.tempBias
mchTropBias=mch+args.tropBias
if args.tempCorr != 0.0:
    tempCorr=mch+args.tempCorr
if args.tropCorr != 0.0:
    tropCorr=mch+args.tropCorr

# Calculate biomass with bias
acdTempBias=0.378*(mchTempBias**2)
acdTempBiasPct=(acdTempBias/acdTemp)*100
acdTropBias=0.844*(mchTropBias**2)
acdTropBiasPct=(acdTropBias/acdTrop)*100
if args.tempCorr != 0.0:
    tempCorrBias=0.378*(tempCorr**2)
    tempCorrBiasPct=(tempCorrBias/acdTemp)*100
if args.tropCorr != 0.0:
    tropCorrBias=0.844*(tropCorr**2)
    tropCorrBiasPct=(tropCorrBias/acdTrop)*100

fig,ax=plt.subplots()
ax.add_patch(Rectangle((5,80),65,40,facecolor='lightgrey'))
ax.plot(mch,acdTempBiasPct,label='Temperate Forest ' + str(args.tempBias) + ' m Bias',color='blue')
ax.plot(mch,acdTropBiasPct,label='Tropical Forest ' + str(args.tropBias) + ' m Bias',color='orange')
if args.tempCorr != 0.0:
    ax.plot(mch,tempCorrBiasPct,label='Temperate Forest ' + str(args.tempCorr) + ' m Bias',linestyle='--',color='blue')
if args.tropCorr != 0.0:
    ax.plot(mch,tropCorrBiasPct,label='Tropical Forest ' + str(args.tropCorr) + ' m Bias',linestyle='--',color='orange')
ax.axhline(y=100,color='black',linestyle='--')
plt.legend()
plt.xlim([5,70])
plt.ylim([0,150])
plt.xlabel('Mean Canopy Height (m)')
plt.ylabel('Biomass (% of Reference Value)')
plt.savefig(args.output + '.png',dpi=300)
plt.close()
plt.clf()
#plt.show()
