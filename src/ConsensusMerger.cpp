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

#include "ConsensusMerger.h"

namespace quandenser {

/* Mcc = mass-charge candidate */
void ConsensusMerger::mergeMccs(std::vector<maracluster::MassChargeCandidate>& allMccs, 
      std::vector<maracluster::MassChargeCandidate>& consensusMccs, int clusterIdx) {
  if (spectrumClusterToConsensusFeatures_.find(clusterIdx) != spectrumClusterToConsensusFeatures_.end()) {
    std::vector<DinosaurFeature>::const_iterator featureIt;
    for (featureIt = spectrumClusterToConsensusFeatures_.find(clusterIdx)->second.begin(); 
          featureIt != spectrumClusterToConsensusFeatures_.find(clusterIdx)->second.end();
          ++featureIt) {
      consensusMccs.push_back(maracluster::MassChargeCandidate(featureIt->charge, 
          featureIt->precMz, maracluster::SpectrumHandler::calcMass(featureIt->precMz, featureIt->charge)));
    }
  }
  /* else {
    maracluster::MSClusterMerge::mergeMccs(allMccs, consensusMccs);
  }*/
}

} /* namespace quandenser */
