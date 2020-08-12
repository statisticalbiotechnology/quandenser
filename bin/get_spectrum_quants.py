'''
Obtain a tab separated file with intensities per run for each consensus 
spectrum, together with the feature-match probabilities from a Quandenser
feature groups output file. Requires Triqler and NumPy to be installed.

Example:

python get_spectrum_quants.py Quandenser_output/Quandenser.feature_groups.tsv --file_list_file file_list.txt

'''

from __future__ import print_function

import os
import sys
import numpy as np
import collections

import triqler.triqler as triqler_base
import triqler.parsers as parsers
import triqler.convert.quandenser as quandenser_convert

def main():
  print('Issued command:', os.path.basename(__file__) + " " + " ".join(map(str, sys.argv[1:])))
  
  args, params = parseArgs()
  
  triqlerInputFile = args.out_file + ".tmp"
  quandenser_convert.convertQuandenserToTriqler(args.file_list_file, args.in_file, [], triqlerInputFile, params)
  
  peptQuantRowMap = collections.defaultdict(list)
  seenSpectra = set()
  runCondPairs = list()
  for i, trqRow in enumerate(parsers.parseTriqlerInputFile(triqlerInputFile)):
    if i % 1000000 == 0:
      print("  Reading row", i)
    
    peptQuantRowMap[trqRow.featureClusterId].append(trqRow)
    if (trqRow.run, trqRow.condition) not in runCondPairs:
      runCondPairs.append((trqRow.run, trqRow.condition))
    
    if not np.isnan(trqRow.searchScore) and trqRow.spectrumId not in seenSpectra:
      seenSpectra.add(trqRow.spectrumId)
  
  params['fileList'], params['groupLabels'], params['groups'] = triqler_base._getFilesAndGroups(runCondPairs)
  
  getPEPFromScore = lambda x : 1.0
  
  _, spectrumQuantRows, intensityDiv = triqler_base._selectBestFeaturesPerRunAndPeptide(
        peptQuantRowMap, getPEPFromScore, params, 
        groupingKey = lambda x : x.spectrumId)
  
  headers = ["consensusSpecId", "featureGroupId"] + list(map(lambda x : "intensity_" + x, params['fileList'])) + list(map(lambda x : "linkPEP_" + x, params['fileList']))
  writer = parsers.getTsvWriter(args.out_file)
  writer.writerow(headers)
  for row in spectrumQuantRows:
    if row.spectrum >= 0 or args.retain_without_spectrum:
      row = [row.spectrum, row.featureGroup] + list(map(lambda x : '%.0f' % x, row.quant)) + list(map(lambda x : '%.5f' % x, row.linkPEP))
      writer.writerow(row)
  
  os.remove(triqlerInputFile)

def parseArgs():
  import argparse
  apars = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  requiredNamed = apars.add_argument_group('required arguments')
  
  apars.add_argument('in_file', default=None, metavar = "IN_FILE",
                     help='''Quandenser output file with feature groups.
                          ''')
  
  requiredNamed.add_argument('--file_list_file', metavar='L', 
                     help='Simple text file with spectrum file names in first column and condition in second column.',
                     required = True)
  
  apars.add_argument('--out_file', default = "spectrum_quants.tsv", metavar='OUT', 
                     help='''Path to triqler input file (writing in TSV format).
                          ''')
  
  apars.add_argument('--min_samples', type=int, default=1, metavar='N', 
                     help='Minimum number of samples a feature group needed to be quantified in.')
                     
  apars.add_argument('--skip_normalization',
                     help='Skip retention-time based intensity normalization.',
                     action='store_true')
  
  apars.add_argument('--retain_without_spectrum',
                     help='Keeps MS1 feature groups without spectrum.',
                     action='store_true')
  
  # ------------------------------------------------
  args = apars.parse_args()
  
  params = dict()
  params['skipNormalization'] = args.skip_normalization
  params['minSamples'] = args.min_samples
  
  params['simpleOutputFormat'] = False
  params["groups"] = []
  params["groupLabels"] = []
  params['retainUnidentified'] = True
  params['plotScatter'] = False
  params['hasLinkPEPs'] = True
  params['writeSpectrumQuants'] = True
  params['decoyPattern'] = 'decoy_'
  
  return args, params
  

if __name__ == "__main__":
   main()
