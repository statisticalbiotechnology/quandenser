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

#include "SpectrumToPrecursorMap.h"

namespace quandenser {

void SpectrumToPrecursorMap::serialize(const std::string& outputFile) { 
  bool append = false;
  std::vector<std::map<int, DinosaurFeatureList> >::const_iterator fileIt;
  int fileIdx;
  for (fileIt = spectrumToPrecursorMap_.begin(), fileIdx = 0; 
       fileIt != spectrumToPrecursorMap_.end(); ++fileIt, ++fileIdx) {
    std::vector<SpectrumFeaturePair> outputVector;
    std::map<int, DinosaurFeatureList>::const_iterator spectrumIt;
    if (Globals::VERB > 1) {
      std::cerr << "Serializing " << fileIt->size() 
        << " spectrum to precursor pairs for file idx " << fileIdx << std::endl;
    }
    for (spectrumIt = fileIt->begin(); spectrumIt != fileIt->end(); ++spectrumIt) {
      DinosaurFeatureList::const_iterator ftIt;
      for (ftIt = spectrumIt->second.begin(); ftIt != spectrumIt->second.end(); ++ftIt) {
        SpectrumFeaturePair specFtPair;
        specFtPair.scanId.fileIdx = fileIdx;
        specFtPair.scanId.scannr = spectrumIt->first;
        specFtPair.feature = *ftIt;
        outputVector.push_back(specFtPair);
      }
    }
    maracluster::BinaryInterface::write(outputVector, outputFile, append);
    append = true;
  }
}

void SpectrumToPrecursorMap::deserialize(const std::string& inputFile) {
  std::vector<SpectrumFeaturePair> inputVector;
  maracluster::BinaryInterface::read(inputFile, inputVector);
  
  std::vector<SpectrumFeaturePair>::iterator pairIt;
  for (pairIt = inputVector.begin(); pairIt != inputVector.end(); ++pairIt) {
    addFeature(pairIt->scanId, pairIt->feature);
  }
}

} /* namespace quandenser */
