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
#include "maracluster/src/BinaryInterface.h"

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

class FeatureToGroupMap {
 public:
  FeatureToGroupMap() : featureIdToGroupId_() {};
  ~FeatureToGroupMap(){};
  
  inline size_t& operator[](const size_t n) { return featureIdToGroupId_[n]; }
  
  inline bool contains(const size_t n) const { 
    return featureIdToGroupId_.find(n) != featureIdToGroupId_.end();
  }
  inline void clear() { featureIdToGroupId_.clear(); }
  
  size_t getGroupId(const size_t n, const std::map<size_t, size_t>& groupIdMap) {
    size_t clusterIdx = featureIdToGroupId_.at(n);
    while (groupIdMap.find(clusterIdx) != groupIdMap.end()) {
      clusterIdx = groupIdMap.at(clusterIdx);
    }
    return clusterIdx;
  }
  
  size_t loadFromFile(const std::string& mapFile) {
    std::vector<std::pair<int, size_t> > featureGroupPairs;
    maracluster::BinaryInterface::read(mapFile, featureGroupPairs);
    for (std::pair<int, size_t>& featureGroupPair : featureGroupPairs) {
      featureIdToGroupId_[featureGroupPair.first] = featureGroupPair.second;
    }
    return featureIdToGroupId_.size();
  }
  
  void saveToFile(const std::string& ftFile, bool append) {
    std::vector<std::pair<int, size_t> > featureGroupPairs(featureIdToGroupId_.begin(), featureIdToGroupId_.end());
    maracluster::BinaryInterface::write(featureGroupPairs, ftFile, append);
  }
 protected:
  std::map<int, size_t> featureIdToGroupId_;
  
};

class FeatureGroups {
 public:
  FeatureGroups(int maxMissingValues, float intensityScoreThreshold) : 
    maxMissingValues_(maxMissingValues), 
    intensityScoreThreshold_(intensityScoreThreshold),
    featureGroups_() {}
  
  void singleLinkClustering(
    const std::vector<std::pair<int, FilePair> >& featureAlignmentQueue,
    std::map<FilePair, std::map<int, FeatureIdxMatch> >& featureMatches,
    const std::string& tmpFilePrefix);
  
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
    FeatureToGroupMap& queryFeatureIdToGroupId,
    FeatureToGroupMap& targetFeatureIdToGroupId,
    std::map<size_t, size_t>& groupIdMap);
  
  void mergeFeatureGroups(const size_t clusterIdx, 
    const size_t mergeInClusterIdx);
  
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
