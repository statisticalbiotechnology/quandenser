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

#ifndef QUANDENSER_SPECTRUMTOPRECURSORMAP_H_
#define QUANDENSER_SPECTRUMTOPRECURSORMAP_H_

#include <vector>
#include <map>

#include "maracluster/src/ScanId.h"
#include "maracluster/src/BinaryInterface.h"

#include "Globals.h"
#include "DinosaurFeature.h"
#include "DinosaurFeatureList.h"

namespace quandenser {

struct SpectrumFeaturePair {
  maracluster::ScanId scanId;
  DinosaurFeature feature;
};

class SpectrumToPrecursorMap {
 public:
  SpectrumToPrecursorMap(const size_t numFiles) : spectrumToPrecursorMap_(numFiles) {}
  
  inline bool isInitialized(maracluster::ScanId& scanId) { 
    return spectrumToPrecursorMap_[scanId.fileIdx][scanId.scannr].isInitialized(); 
  }
  inline void setInitialized(maracluster::ScanId& scanId) { 
    spectrumToPrecursorMap_[scanId.fileIdx][scanId.scannr].setInitialized(); 
  }
  inline void addFeature(maracluster::ScanId& scanId, DinosaurFeature ft) { 
    spectrumToPrecursorMap_[scanId.fileIdx][scanId.scannr].push_back(ft); 
  }
  
  inline std::vector<DinosaurFeature>::const_iterator getBegin(maracluster::ScanId& scanId) { 
    return spectrumToPrecursorMap_[scanId.fileIdx][scanId.scannr].begin();
  }
  inline std::vector<DinosaurFeature>::const_iterator getEnd(maracluster::ScanId& scanId) {
    return spectrumToPrecursorMap_[scanId.fileIdx][scanId.scannr].end();
  }
  
  void serialize(const std::string& outputFile);
  void deserialize(const std::string& inputFile);
  
 protected:
  std::vector<std::map<int, DinosaurFeatureList> > spectrumToPrecursorMap_;
};

} /* namespace quandenser */

#endif /* QUANDENSER_SPECTRUMTOPRECURSORMAP_H_ */
