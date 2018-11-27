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

#ifndef QUANDENSER_FEATUREGROUPS_H_
#define QUANDENSER_FEATUREGROUPS_H_

#include <vector>
#include <map>
#include <set>
#include <cmath>

#include "maracluster/src/SpectrumFileList.h"
#include "maracluster/src/ScanId.h"

#include "DinosaurFeatureList.h"
#include "SimilarityMatrix.h"
#include "FeatureAlignment.h"

namespace quandenser {

typedef maracluster::ScanId FeatureId;

struct FeatureIdMatch {
  FeatureIdMatch() : queryFeatureId(-1, -1), targetFeatureId(-1, -1), 
                     posteriorErrorProb(1.0) {}
  FeatureIdMatch(FeatureId qFtId, FeatureId tFtId, float pep) : 
    queryFeatureId(qFtId), targetFeatureId(tFtId), 
    posteriorErrorProb(pep) {}
  
  FeatureId queryFeatureId, targetFeatureId;
  float posteriorErrorProb;
};

bool operator<(const FeatureIdMatch& l, const FeatureIdMatch& r);

class FeatureGroups {
 public:
  FeatureGroups(int maxMissingValues, float intensityScoreThreshold) : 
    maxMissingValues_(maxMissingValues), 
    intensityScoreThreshold_(intensityScoreThreshold),
    featureGroups_() {}
  
  void singleLinkClustering(
    const std::vector<std::pair<int, FilePair> >& featureAlignmentQueue,
    std::map<FilePair, std::map<int, FeatureIdxMatch> >& featureMatches);
  
  void filterConsensusFeatures(
    const std::vector<DinosaurFeatureList>& allFeatures,
    std::map<FeatureId, std::vector<int> >& featureToSpectrumCluster);
  
  void printFeatureGroups(
    const std::string& featureGroupsOutFile,
    std::vector<DinosaurFeatureList>& allFeatures,
    const std::map<FeatureId, std::vector<int> >& featureToSpectrumCluster,
    std::map<int, std::vector<DinosaurFeature> >& spectrumClusterToConsensusFeatures);
    
 protected:  
  int maxMissingValues_;
  float intensityScoreThreshold_;
  
  std::vector<std::vector<FeatureIdMatch> > featureGroups_;
  
  void addToFeatureGroup(const FeatureId& queryFeatureId,
    const FeatureId& targetFeatureId, float posteriorErrorProb,
    std::map<FeatureId, size_t>& featureIdToGroupId);
  
  void addConsensusFeatureToSpectrumClusters(
    const DinosaurFeature& consensusFeature,
    std::map<FeatureId, std::map<int, float> >& spectrumClusterLinkPEPs,
    std::map<int, std::vector<DinosaurFeature> >& spectrumClusterToConsensusFeatures,
    std::map<int, size_t>& spectrumClusterIdxOffsets);
  
  void printFeatureGroup(
    const std::vector<DinosaurFeature>& features,
    std::map<FeatureId, std::map<int, float> >& spectrumClusterLinkPEPs,
    std::map<int, size_t>& spectrumClusterIdxOffsets,
    std::ostream& dataStream);
  
  void propagateSpectrumClusterIds(
    const std::set<FeatureId>& featureIds, 
    const std::map<FeatureId, std::vector<int> >& featureToSpectrumCluster, 
    SimilarityMatrix<FeatureId>& simMatrix,
    std::map<FeatureId, std::map<int, float> >& spectrumClusterLinkPEPs);
  
  void updateLinkPEPs(const FeatureId rootFeatureId,
    const std::vector<int>& clusterIdxs,
    const std::map<FeatureId, float>& similarities,
    std::map<FeatureId, std::map<int, float> >& spectrumClusterLinkPEPs);
  
};

} /* namespace quandenser */

#endif /* QUANDENSER_FEATUREGROUPS_H_ */
