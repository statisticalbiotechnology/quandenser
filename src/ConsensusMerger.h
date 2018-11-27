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

#ifndef QUANDENSER_CONSENSUSMERGER_H_
#define QUANDENSER_CONSENSUSMERGER_H_

#include <vector>
#include <map>
#include <sstream>

#include "maracluster/src/MSFileMerger.h"
#include "maracluster/src/MassChargeCandidate.h"

#include "Globals.h"
#include "DinosaurFeature.h"

namespace quandenser {

class ConsensusMerger : public maracluster::MSFileMerger {
 public:
  ConsensusMerger(std::string& spectrumOutFN, 
      const std::map<int, std::vector<DinosaurFeature> >& map) : 
    maracluster::MSFileMerger(spectrumOutFN), spectrumClusterToConsensusFeatures_(map) {}
  
 protected:
  /* Mcc = mass-charge candidate */
  void mergeMccs(std::vector<maracluster::MassChargeCandidate>& allMccs, 
      std::vector<maracluster::MassChargeCandidate>& consensusMccs, int clusterIdx);
  
  const std::map<int, std::vector<DinosaurFeature> >& spectrumClusterToConsensusFeatures_;
    
};

} /* namespace quandenser */

#endif /* QUANDENSER_CONSENSUSMERGER_H_ */
