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

#ifndef QUANDENSER_MARACLUSTERADAPTER_H_
#define QUANDENSER_MARACLUSTERADAPTER_H_

#include <vector>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "maracluster/src/MaRaCluster.h"

#include "SpectrumFiles.h"
#include "DinosaurFeature.h"
#include "DinosaurFeatureList.h"
#include "SpectrumToPrecursorMap.h"
#include "ConsensusMerger.h"

namespace quandenser {

class MaRaClusterAdapter : public maracluster::MaRaCluster {
 public:
  MaRaClusterAdapter(const std::vector<DinosaurFeatureList>& featureLists,
      SpectrumToPrecursorMap& specToPrecMap, std::string spectrumOutputFile) :
        MaRaCluster(), featureLists_(featureLists), 
        spectrumToPrecursorMap_(specToPrecMap) {
    spectrumOutFN_ = spectrumOutputFile;
  }
  
  bool parseOptions(const std::vector<std::string>& maraclusterArgs) {
    std::vector<const char*> maraclusterArgv;
    for (std::vector<std::string>::const_iterator it = maraclusterArgs.begin();
       it != maraclusterArgs.end(); ++it) {
      maraclusterArgv.push_back(it->c_str());
    }
    return maracluster::MaRaCluster::parseOptions(
        maraclusterArgs.size(), 
        (char**)&maraclusterArgv.front());
  }
  
  int mergeSpectra(); // overrides base function
  
  inline std::map<int, std::vector<DinosaurFeature> >& getSpectrumClusterToConsensusFeatures() { 
    return spectrumClusterToConsensusFeatures_; 
  }
  inline const std::map<int, std::vector<DinosaurFeature> >& getSpectrumClusterToConsensusFeatures() const { 
    return spectrumClusterToConsensusFeatures_; 
  }
  
  std::string getClusterFileName() {
    int negIntThreshold = -1*static_cast<int>(clusterThresholds_.back());
    return outputFolder_ + "/" + fnPrefix_ + ".clusters_p" + 
               boost::lexical_cast<std::string>(negIntThreshold) + ".tsv"; 
  }
  
 protected:
  int createIndex(); // overrides base function
  
  const std::vector<DinosaurFeatureList>& featureLists_;
  SpectrumToPrecursorMap& spectrumToPrecursorMap_;
  
  std::map<int, std::vector<DinosaurFeature> > spectrumClusterToConsensusFeatures_;
    
};

} /* namespace quandenser */

#endif /* QUANDENSER_MARACLUSTERADAPTER_H_ */
