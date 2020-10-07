
# Python3 compatibility
from __future__ import print_function, division, absolute_import

# Load DmpSoftware and ROOT
from ROOT import gSystem
gSystem.Load('libDmpEvent.so')
from ROOT import DmpChain, DmpEvent

# Python libraries
import matplotlib   # Plotting library
matplotlib.use('Agg')  # Non-interactive mode
import matplotlib.pyplot as plt

# Very high energy proton file

inputFile1 = '/beegfs/dampe/prod/MC/reco/v6r0p0/allProton-v6r0p0_10TeV_100TeV-HP-p6/allProton-v6r0p0_10TeV_100TeV-HP-p6.noOrb.607416.reco.root'
inputFile2 = '/beegfs/dampe/prod/MC/reco/v6r0p0/allProton-v6r0p0_10TeV_100TeV-HP-p6/allProton-v6r0p0_10TeV_100TeV-HP-p6.noOrb.601082.reco.root'

# DmpChain object: http://dpnc.unige.ch/dampe/doxygen/DmpSoftware-6-0-0/classDmpChain.html
# DmpChain is inherited from standard ROOT TChain
dc = DmpChain("CollectionTree") 
dc.Add(inputFile1)
dc.Add(inputFile2)  # Can chain several files together

nevents = dc.GetEntries()

storageArray = []  # Empty list

for i in range(nevents):

	pev = dc.GetDmpEvent(i)
	# DmpEvent object: http://dpnc.unige.ch/dampe/doxygen/DmpSoftware-6-0-0/classDmpEvent.html
	
	BGO_calorimeter = pev.pEvtBgoRec()
	measuredEnergy = BGO_calorimeter.GetTotalEnergy()

	storageArray.append(measuredEnergy)
	# We could of course write directly:
	# storageArray.append(pev.pEvtBgoRec().GetTotalEnergy() )


# Plot the energy histogram
fig1 = plt.figure()
plt.hist(storageArray,50,histtype='step')
plt.xscale('log')  # Limitation: the binning is not logarithmic
plt.yscale('log')
plt.xlabel('Reconstructed energy')
plt.ylabel('Event count')
plt.title('allProton-v6r0p0_10TeV_100TeV-HP-p6')
plt.savefig('histogram.png')
plt.close(fig1)
