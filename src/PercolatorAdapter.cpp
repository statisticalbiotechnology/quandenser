/******************************************************************************  
  Copyright 2015-2017 Matthew The <matthew.the@scilifelab.se>
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
 ******************************************************************************/

#include "PercolatorAdapter.h"

namespace quandenser {

PercolatorAdapter::PercolatorAdapter() : Caller(), setHandler_(0u) {}
  
void PercolatorAdapter::process() {
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  if (Globals::VERB > 0) {
    cerr << extendedGreeter(startTime);
  }
  
  if (Globals::VERB > 2) {
    std::cerr << "FeatureNames::getNumFeatures(): "<< FeatureNames::getNumFeatures() << endl;
  }
  
  setHandler_.normalizeFeatures(pNorm_);
  
  // Copy feature data pointers to Scores object
  Scores allScores(useMixMax_);
  allScores.populateWithPSMs(setHandler_);
  
  CrossValidation crossValidation(quickValidation_, reportEachIteration_, 
                                  testFdr_, selectionFdr_, initialSelectionFdr_, selectedCpos_, 
                                  selectedCneg_, numIterations_, useMixMax_,
                                  nestedXvalBins_, trainBestPositive_,
                                  numThreads_, skipNormalizeScores_);
  int firstNumberOfPositives = crossValidation.preIterationSetup(allScores, pCheck_, pNorm_, setHandler_.getFeaturePool());
  if (Globals::VERB > 0) {
    cerr << "Found " << firstNumberOfPositives << " test set positives with q<"
        << testFdr_ << " in initial direction" << endl;
  }
  
  time_t procStart;
  clock_t procStartClock = clock();
  time(&procStart);
  double diff = difftime(procStart, startTime);
  if (Globals::VERB > 1) cerr << "Reading in data and feature calculation took "
      << ((double)(procStartClock - startClock)) / (double)CLOCKS_PER_SEC
      << " cpu seconds or " << diff << " seconds wall clock time." << endl;
  
  if (tabOutputFN_.length() > 0) {
    setHandler_.writeTab(tabOutputFN_, pCheck_);
  }
  
  // Do the SVM training
  crossValidation.train(pNorm_);
  
  if (weightOutputFN_.size() > 0) {
    ofstream weightStream(weightOutputFN_.c_str(), ios::out);
    crossValidation.printAllWeights(weightStream, pNorm_);
    weightStream.close();
  }
  
  // Calculate the final SVM scores and clean up structures
  crossValidation.postIterationProcessing(allScores, pCheck_);
  
  // calculate psms level probabilities TDA or TDC
  bool isUniquePeptideRun = false;
  calculatePSMProb(allScores, isUniquePeptideRun, procStart, procStartClock, diff);
  
  /* These are all example of why you should not use the singleton pattern... */
  DataSet::resetFeatureNames();
  Normalizer::resetNormalizer();
}
  
void PercolatorAdapter::init(const std::vector<std::pair<std::string, double> >& percolatorFeatures) {
  // read feature names and initial values that are present in feature descriptions
  FeatureNames& featureNames = DataSet::getFeatureNames();
  std::vector<std::pair<std::string, double> >::const_iterator featureIt;
  std::vector<double> init_values(percolatorFeatures.size());
  bool hasDefaultValues = false;
  int i = 0;
  for (featureIt = percolatorFeatures.begin(); featureIt != percolatorFeatures.end(); ++featureIt, ++i) {
    featureNames.insertFeature(featureIt->first);
    if (featureIt->second != 0.0) {
      hasDefaultValues = true;
      init_values[i] = featureIt->second;
      if (Globals::VERB > 3) {
        std::cerr << "Initial direction for " << featureIt->first << " is " << 
                     featureIt->second << std::endl;
      }
    }
  }
  bool calcDoc = false;
  featureNames.initFeatures(calcDoc);
  
  setHandler_.getFeaturePool().createPool(FeatureNames::getNumFeatures());
  
  DataSet* targetSet = new DataSet();
  assert(targetSet);
  targetSet->setLabel(1);
  DataSet* decoySet = new DataSet();
  assert(decoySet);
  decoySet->setLabel(-1);    

  setHandler_.push_back_dataset(targetSet);
  setHandler_.push_back_dataset(decoySet);
  
  pCheck_ = SanityCheck::initialize("");
  assert(pCheck_);
  pCheck_->checkAndSetDefaultDir();
  if (hasDefaultValues) pCheck_->addDefaultWeights(init_values);
  bool concatenatedSearch = true;
  pCheck_->setConcatenatedSearch(concatenatedSearch);
}

void PercolatorAdapter::addPsm(const std::string& specId, int label, int scannr, 
    const std::string& peptide, const std::vector<double>& features) {
  PSMDescription* myPsm = new PSMDescription();
  
  // adding n-term and c-term residues to peptide
  myPsm->peptide = peptide;

  myPsm->setId(specId);
  myPsm->scan = scannr;
  myPsm->expMass = 0.0;
  myPsm->calcMass = 0.0;

  myPsm->features = setHandler_.getFeaturePool().allocate();
  
  std::vector<double>::const_iterator ftIt;
  int i = 0;
  for (ftIt = features.begin(); ftIt != features.end(); ++ftIt, ++i) {
    myPsm->features[i] = *ftIt;
  }
  
  setHandler_.getSubsetFromLabel(label)->registerPsm(myPsm);
}

} /* namespace quandenser */
