from __future__ import print_function

import csv
import sys
import os
import numpy as np
import collections
import re
import getopt

from triqler import parsers
from triqler import triqler

import normalize_intensities as normalize
import percolator

def main(argv):
  helpMessage = 'prepare_input.py -l <batch_file_list> -f <cluster_quant_file> -i <pout_file_target,pout_file_decoy> -q <triqler_input_file>'
  
  simpleOutputFormat = False
  try:
    opts, args = getopt.getopt(argv,"hi:l:o:f:q:e:s",["ifile=","filelist=","featurefile=","simpleoutput="])
  except getopt.GetoptError:
    print(helpMessage)
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print(helpMessage)
      sys.exit()
    elif opt in ("-i", "--ifile"):
      psmsOutputFiles = arg.split(",")
    elif opt in ("-l", "--filelist"):
      fileListFile = arg
    elif opt in ("-f", "--featurefile"):
      clusterQuantFile = arg
    elif opt in ("-s", "--simpleoutput"):
      simpleOutputFormat = True
    elif opt in ("-q"):
      peptQuantRowFile = arg

  clusterQuantFileNormalized = clusterQuantFile.replace(".tsv", ".normalized.tsv")
  if not os.path.isfile(clusterQuantFileNormalized):
    normalize.normalizeIntensitiesRtimeBased(clusterQuantFile, clusterQuantFileNormalized)
  clusterQuantFile = clusterQuantFileNormalized
  
  params = dict()
  params['simpleOutputFormat'] = simpleOutputFormat
  fileList, params['groups'], params['groupLabels'] = parsers.parseFileList(fileListFile)
  specToPeptideMap = parsePsmsPoutFiles(psmsOutputFiles)
  
  printTriqlerInputFile(fileList, clusterQuantFile, peptQuantRowFile, specToPeptideMap, params)
    
def parsePsmsPoutFiles(psmsOutputFiles):
  specToPeptideMap = collections.defaultdict(list)
  for psmsOutputFile in psmsOutputFiles:
    for psm in percolator.parsePsmsPout(psmsOutputFile):
      specToPeptideMap[psm.scannr] = (psm.peptide, psm.PEP, psm.proteins, psm.svm_score, psm.charge)
  return lambda spectrumIdx : specToPeptideMap.get(spectrumIdx, getDefaultPeptideHit(spectrumIdx))

def getDefaultPeptideHit(spectrumIdx):
  return ("NA", 1.0, ["NA"], np.nan, -1) # psm.peptide, psm.PEP, psm.proteins, psm.svm_score, psm.charge
  
def parsePeptideLinkPEP(peptLinkPEP):
  spectrumIdx, linkPEP = peptLinkPEP.split(";")
  return int(spectrumIdx), float(linkPEP)

def printTriqlerInputFile(fileList, clusterQuantFile, quantRowFile, specToPeptideMap, params):
  print("Parsing cluster quant file")
  
  fileNameConditionPairs = [[x.split("/")[-1], parsers.getGroupLabel(idx, params['groups'], params['groupLabels'])] for idx, x in enumerate(fileList)]
  
  writer = csv.writer(open(quantRowFile, 'w'), delimiter = '\t')
  if params['simpleOutputFormat']:
    writer.writerow(parsers.TriqlerSimpleInputRowHeaders)
  else:
    writer.writerow(parsers.TriqlerInputRowHeaders)
  
  featureClusterRows = list()
  spectrumToFeatureMatch = dict() # stores the best peptideQuantRow per (peptide, spectrumIdx)-pair
  for featureClusterIdx, featureCluster in enumerate(parsers.parseFeatureClustersFile(clusterQuantFile)):
    if featureClusterIdx % 10000 == 0:
      print("Processing feature group", featureClusterIdx + 1)
    
    rows = list()
    for pc in featureCluster:
      fileIdx = int(pc.fileName)
      for peptLinkPEP in pc.peptLinkPEPs.split(","):
        spectrumIdx, linkPEP = parsePeptideLinkPEP(peptLinkPEP)
        peptide, identPEP, proteins, searchScore, charge = specToPeptideMap(spectrumIdx)
        if pc.intensity > 0.0 and linkPEP < 1.0:
          # run condition charge spectrumId linkPEP featureClusterId search_score intensity peptide proteins
          run, condition = fileNameConditionPairs[fileIdx]
          row = parsers.TriqlerInputRow(run, condition, charge, spectrumIdx, linkPEP, featureClusterIdx, searchScore, pc.intensity, peptide, proteins)
          rows.append(row)
    
    newRows = list()
    rows = sorted(rows, key = lambda x : (x.run, x.spectrumId, x.linkPEP, -1*x.searchScore))
    prevKey = (-1, -1)
    bestSearchScore = -1e9
    for row in rows:
      if prevKey == (row.run, row.spectrumId):
        if row.searchScore > bestSearchScore:
          bestSearchScore = row.searchScore
          newRows.append(row)
      else:
        newRows.append(row)
        prevKey = (row.run, row.spectrumId)
    
    for row in newRows:
      if params['simpleOutputFormat']:
        writer.writerow(row.toSimpleList())
      else:
        writer.writerow(row.toList())
  
def combinePEPs(linkPEP, identPEP):
  return 1.0 - (1.0 - linkPEP)*(1.0 - identPEP)

if __name__ == "__main__":
   main(sys.argv[1:])
