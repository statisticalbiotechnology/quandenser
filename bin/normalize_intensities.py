from __future__ import print_function

import sys
import csv
import bisect
from collections import defaultdict

import numpy as np

from triqler import parsers

def normalizeIntensitiesRtimeBased(clusterQuantExtraFile, clusterQuantExtraNormalizedFile, plotScatter = False, plotRunningAverage = False):
  intensities = getIntensityFactorPairs(clusterQuantExtraFile)
  
  if plotScatter:
    plotFactorScatter(intensities)
  
  rTimeFactorArrays = getFactorArrays(intensities)
  
  if plotRunningAverage:
    plotFactorRunningAverage(rTimeFactorArrays)
  
  normalizeIntensitiesWithFactorArrays(clusterQuantExtraFile, rTimeFactorArrays, clusterQuantExtraNormalizedFile)
  
def getIntensityFactorPairs(clusterQuantExtraFile, minRunsObservedIn = 4):
  precClusters = parsers.parseFeatureClustersFile(clusterQuantExtraFile)
  intensities = defaultdict(list)
  for i, precCluster in enumerate(precClusters):
    intensity, rtime = dict(), dict()
    sortedPrecursors = sorted(precCluster, key = lambda x : (x.fileName, -1.0 * x.intensity))
    for row in sortedPrecursors:
      if not np.isnan(row.intensity) and not "_decoy" in row.fileName:
        if row.fileName not in intensity: # and row.charge == 2:
          intensity[row.fileName] = (np.log2(row.intensity), row.rTime)
    if len(intensity) >= minRunsObservedIn:
      keys = sorted(intensity.keys())
      masterKey = keys[0]
      factor0 = sum([intensity[masterKey][0] - x[0] for x in intensity.values()]) / len(keys)
      intensities[masterKey].append((intensity[masterKey][1], factor0))
      for key in keys[1:]:
        factori = intensities[masterKey][-1][1] - (intensity[masterKey][0] - intensity[key][0])
        intensities[key].append((intensity[key][1], factori))
    if (i+1) % 10000 == 0:
      print("Processing cluster", i+1)
  return intensities

# returns running averages of factors
def getFactorArrays(intensities, N = 2000):
  factorArrays = defaultdict(list)
  for i, key in enumerate(sorted(intensities.keys())):
    print("Calculating factor array", i+1)
    intensities[key] = sorted(intensities[key], key = lambda x : x[0])
    rTimes = [x[0] for x in intensities[key]]
    factors = [x[1] for x in intensities[key]]
    factorArrays[key] = zip(rTimes[int(N/2):int(-N/2)], runningMean(factors, N))
  return factorArrays

def runningMean(x, N):
  cumsum = np.cumsum(np.insert(x, 0, 0)) 
  return (cumsum[N:] - cumsum[:-N]) / N 
  
def normalizeIntensitiesWithFactorArrays(clusterQuantExtraFile, rTimeFactorArrays, clusterQuantExtraNormalizedFile):
  rTimeArrays, factorArrays = dict(), dict()
  for key in rTimeFactorArrays:
    rTimeArrays[key], factorArrays[key] = zip(*rTimeFactorArrays[key])
  print("Writing", clusterQuantExtraNormalizedFile)
  writer = csv.writer(open(clusterQuantExtraNormalizedFile, 'w'), delimiter = '\t')
  precClusters = parsers.parseFeatureClustersFile(clusterQuantExtraFile)
  for i, precCluster in enumerate(precClusters):
    for row in precCluster:
      rTimeIndex = min([bisect.bisect_left(rTimeArrays[row.fileName], row.rTime), len(rTimeArrays[row.fileName]) - 1])
      outRow = list(row[1:])
      outRow[4] = row.intensity / (2 ** factorArrays[row.fileName][rTimeIndex])
      writer.writerow(outRow)
    writer.writerow([])
    if (i+1) % 10000 == 0:
      print("Writing cluster", i+1)
  
def plotFactorScatter(intensities):
  import scatter
  import matplotlib.pyplot as plt
  scatter.prepareSubplots()
  for i, key in enumerate(sorted(intensities.keys())):
    plt.subplot((len(intensities)-1) / 4 + 1,4,i+1)
    print("Calculating density", i+1)
    rTimes = [x[0] for x in intensities[key] if abs(x[1]) < 2]
    factors = [x[1] for x in intensities[key] if abs(x[1]) < 2]
    scatter.plotDensityScatterHist(rTimes, factors, bins = [20,20])
    plt.plot([min(rTimes),max(rTimes)],[0,0],'k:')
    plt.title(key.replace("JD_06232014_","").replace("_","\_"))
    plt.xlabel("Retention time")
    plt.ylabel("log2(int0/int1)")
    plt.ylim([-2,2])
  plt.show()

def plotFactorRunningAverage(factorArrays):
  import scatter
  import matplotlib.pyplot as plt
  scatter.prepareSubplots()
  for i, key in enumerate(sorted(factorArrays.keys())):
    plt.subplot((len(factorArrays)-1) / 4 + 1,4,i+1)
    plt.title(key.replace("JD_06232014_","").replace("_","\_"))
    plt.xlabel("Retention time")
    plt.ylabel("log2(int0/int1)")
    plt.ylim([-2,2])
    plt.plot(*zip(*factorArrays[key]))
  plt.show()

if __name__ == "__main__":
   main(sys.argv[1:])
