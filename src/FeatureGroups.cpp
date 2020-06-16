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

#include "FeatureGroups.h"
#include "boost/lexical_cast.hpp"

namespace quandenser {

bool operator<(const FeatureIdMatch& l, const FeatureIdMatch& r) {
  return (l.posteriorErrorProb < r.posteriorErrorProb) 
    || (l.posteriorErrorProb == r.posteriorErrorProb &&   
        l.queryFeatureId < r.queryFeatureId) 
    || (l.posteriorErrorProb == r.posteriorErrorProb && 
        l.queryFeatureId == r.queryFeatureId && 
        l.targetFeatureId < r.targetFeatureId);
}

void FeatureGroups::singleLinkClustering(
    const std::vector<std::pair<int, FilePair> >& featureAlignmentQueue,
    std::map<FilePair, std::map<int, FeatureIdxMatch> >& featureMatches,
    const std::string& tmpFilePrefixGroup,
    const std::string& tmpFilePrefixAlign) {
  size_t numFiles = featureAlignmentQueue.size() / 2 + 1;
  std::vector<FeatureToGroupMap> featureIdToGroupId(numFiles);
  std::map<size_t, size_t> groupIdMap;
  
  if (!tmpFilePrefixGroup.empty()) {
    for (size_t i = 0; i < numFiles; ++i) {
      std::string fileName = tmpFilePrefixGroup + boost::lexical_cast<std::string>(i) + ".dat";
      remove(fileName.c_str());
    }
  }
  
  std::vector<std::pair<int, FilePair> >::const_iterator filePairIt;
  int filePairCnt = 0;
  for (filePairIt = featureAlignmentQueue.begin(); 
       filePairIt != featureAlignmentQueue.end(); 
       ++filePairIt, ++filePairCnt) {
    FilePair filePair = filePairIt->second;
    int fileIdx1 = filePair.fileIdx1;
    int fileIdx2 = filePair.fileIdx2;
    
    FeatureToGroupMap& queryFeatureIdToGroupId = featureIdToGroupId.at(fileIdx1);
    FeatureToGroupMap& targetFeatureIdToGroupId = featureIdToGroupId.at(fileIdx2);
    
    if (!tmpFilePrefixGroup.empty()) {
      queryFeatureIdToGroupId.loadFromFile(tmpFilePrefixGroup + boost::lexical_cast<std::string>(fileIdx1) + ".dat");
      targetFeatureIdToGroupId.loadFromFile(tmpFilePrefixGroup + boost::lexical_cast<std::string>(fileIdx2) + ".dat");
    }
    
    if (!tmpFilePrefixAlign.empty()) {
      std::string matchesFileName = tmpFilePrefixAlign + "/matches."  + 
        boost::lexical_cast<std::string>(filePair.fileIdx1) + "." + 
        boost::lexical_cast<std::string>(filePair.fileIdx2) + ".dat";
      FeatureAlignment::loadFromFile(matchesFileName, featureMatches[filePair]);
    }
    
    if (Globals::VERB > 2) {
      std::cerr << "Processing pair " << filePairCnt+1 << "/" 
          << featureAlignmentQueue.size() 
          << ": " << fileIdx1 << "->" << fileIdx2 << " (" 
          << featureMatches[filePair].size() << " links)" << std::endl;
    }
    
    std::map<int, FeatureIdxMatch>::const_iterator ftIdxMatchIt;
    for (ftIdxMatchIt = featureMatches[filePair].begin(); 
         ftIdxMatchIt != featureMatches[filePair].end();
         ++ftIdxMatchIt) {
      int queryFeatureIdx = ftIdxMatchIt->first;
      int targetFeatureIdx = ftIdxMatchIt->second.targetFeatureIdx;
      float posteriorErrorProb = ftIdxMatchIt->second.posteriorErrorProb;
      
      /* placeholders have a PEP of 1.0 and will destroy the score, deal with 
         this some other way downstream... */
      if (posteriorErrorProb == 1.0) posteriorErrorProb = 0.0;
      
      FeatureId queryFeatureId(fileIdx1, queryFeatureIdx);
      FeatureId targetFeatureId(fileIdx2, targetFeatureIdx);
      
      addToFeatureGroup(queryFeatureId, targetFeatureId, posteriorErrorProb,
          queryFeatureIdToGroupId, targetFeatureIdToGroupId, groupIdMap);
    }
    
    if (!tmpFilePrefixGroup.empty()) {
      bool append = false;
      queryFeatureIdToGroupId.saveToFile(tmpFilePrefixGroup + boost::lexical_cast<std::string>(fileIdx1) + ".dat", append);
      targetFeatureIdToGroupId.saveToFile(tmpFilePrefixGroup + boost::lexical_cast<std::string>(fileIdx2) + ".dat", append);
      
      queryFeatureIdToGroupId.clear();
      targetFeatureIdToGroupId.clear();
    }
    
    featureMatches[filePair].clear();
  }
  
  if (!tmpFilePrefixGroup.empty()) {
    for (size_t i = 0; i < numFiles; ++i) {
      std::string fileName = tmpFilePrefixGroup + boost::lexical_cast<std::string>(i) + ".dat";
      remove(fileName.c_str());
    }
  }
}

void FeatureGroups::addToFeatureGroup(const FeatureId& queryFeatureId,
    const FeatureId& targetFeatureId, float posteriorErrorProb,
    FeatureToGroupMap& queryFeatureIdToGroupId,
    FeatureToGroupMap& targetFeatureIdToGroupId,
    std::map<size_t, size_t>& groupIdMap) {
  size_t clusterIdx = 0u;
  int queryFeatureIdx = queryFeatureId.scannr;
  int targetFeatureIdx = targetFeatureId.scannr;
  
  if (queryFeatureIdToGroupId.contains(queryFeatureIdx)
       && targetFeatureIdToGroupId.contains(targetFeatureIdx)) {
    clusterIdx = queryFeatureIdToGroupId.getGroupId(queryFeatureIdx, groupIdMap);
    size_t mergeInClusterIdx = targetFeatureIdToGroupId.getGroupId(targetFeatureIdx, groupIdMap);
    
    if (mergeInClusterIdx != clusterIdx) {
      groupIdMap[mergeInClusterIdx] = clusterIdx;
      mergeFeatureGroups(clusterIdx, mergeInClusterIdx);
    }
  } else if (queryFeatureIdToGroupId.contains(queryFeatureIdx)) {
    clusterIdx = queryFeatureIdToGroupId.getGroupId(queryFeatureIdx, groupIdMap);
  } else if (targetFeatureIdToGroupId.contains(targetFeatureIdx)) {
    clusterIdx = targetFeatureIdToGroupId.getGroupId(targetFeatureIdx, groupIdMap);
  } else {
    clusterIdx = featureGroups_.size();
    featureGroups_.push_back(std::vector<FeatureIdMatch>());
  }
  
  queryFeatureIdToGroupId[queryFeatureIdx] = clusterIdx;
  targetFeatureIdToGroupId[targetFeatureIdx] = clusterIdx;
  
  featureGroups_[clusterIdx].push_back(FeatureIdMatch(
      queryFeatureId, targetFeatureId, posteriorErrorProb));
}

void FeatureGroups::mergeFeatureGroups(const size_t clusterIdx, 
    const size_t mergeInClusterIdx) {
  featureGroups_[clusterIdx].insert(
      featureGroups_[clusterIdx].end(), 
      featureGroups_[mergeInClusterIdx].begin(), 
      featureGroups_[mergeInClusterIdx].end() );
  featureGroups_[mergeInClusterIdx].clear();
}

void FeatureGroups::filterConsensusFeatures(
    const std::vector<DinosaurFeatureList>& allFeatures,
    std::map<FeatureId, std::vector<int> >& featureToSpectrumCluster) {
  std::map<int, std::vector< std::pair<float, std::vector<FeatureId> > > > intensityScores;
  
  std::vector<std::vector<FeatureIdMatch> >::iterator ftGroupIt;
  for (ftGroupIt = featureGroups_.begin(); ftGroupIt != featureGroups_.end(); ++ftGroupIt) {
    if (ftGroupIt->empty()) continue;
    
    std::set<FeatureId> featureIds;
    
    std::vector<FeatureIdMatch>::const_iterator ftMatchIt;
    for (ftMatchIt = ftGroupIt->begin(); ftMatchIt != ftGroupIt->end(); ++ftMatchIt) {
      featureIds.insert(ftMatchIt->queryFeatureId);
      featureIds.insert(ftMatchIt->targetFeatureId);
    }
    
    size_t numRuns = allFeatures.size();
    std::map<int, std::pair<std::vector<float>, std::vector<FeatureId> > > maxIntensity;
    for (std::set<FeatureId>::const_iterator ftIt = featureIds.begin(); 
         ftIt != featureIds.end(); ++ftIt) {
      if (featureToSpectrumCluster.find(*ftIt) != featureToSpectrumCluster.end() &&
          !featureToSpectrumCluster.find(*ftIt)->second.empty()) {
        
        int fileIdx = ftIt->fileIdx;
        int featureIdx = ftIt->scannr;
        
        const DinosaurFeature& ft = allFeatures.at(fileIdx).at(featureIdx);
        
        std::vector<int> clusterIdxs(featureToSpectrumCluster.find(*ftIt)->second);
        std::vector<int>::const_iterator clusterIdxIt;
        for (clusterIdxIt = clusterIdxs.begin(); 
             clusterIdxIt != clusterIdxs.end(); ++clusterIdxIt) {
          if (maxIntensity.find(*clusterIdxIt) == maxIntensity.end()) {
            maxIntensity[*clusterIdxIt].first = std::vector<float>(numRuns);
          }
          maxIntensity[*clusterIdxIt].second.push_back(*ftIt);
          if (ft.intensity > maxIntensity[*clusterIdxIt].first[fileIdx]) {
            maxIntensity[*clusterIdxIt].first[fileIdx] = ft.intensity;
          }
        }
      }
    }
    
    std::map<int, std::pair<std::vector<float>, std::vector<FeatureId> > >::const_iterator clusterIt;
    for (clusterIt = maxIntensity.begin(); clusterIt != maxIntensity.end(); ++clusterIt) {
      float intensityScore = 0.0;
      std::vector<float>::const_iterator intensityIt;
      for (intensityIt = clusterIt->second.first.begin(); 
           intensityIt != clusterIt->second.first.end(); ++intensityIt) {
        /* TODO: what to do with intensities between (0,1) which will give negative log scores? */
        if (*intensityIt > 0.0) {
          intensityScore += log(*intensityIt);
        }
      }
      intensityScores[clusterIt->first].push_back(std::make_pair(intensityScore, clusterIt->second.second));
    }
  }
  
  featureToSpectrumCluster.clear();
  std::map<int, std::vector< std::pair<float, std::vector<FeatureId> > > >::iterator clusterIt;
  for (clusterIt = intensityScores.begin(); clusterIt != intensityScores.end(); ++clusterIt) {
    int clusterIdx = clusterIt->first;
    std::sort(clusterIt->second.rbegin(), clusterIt->second.rend());
    
    float maxScore = clusterIt->second.begin()->first;
    std::vector< std::pair<float, std::vector<FeatureId> > >::const_iterator ftGroupIt;
    for (ftGroupIt = clusterIt->second.begin(); ftGroupIt != clusterIt->second.end(); ++ftGroupIt) {
      if (ftGroupIt->first >= maxScore * intensityScoreThreshold_) {
        std::vector<FeatureId>::const_iterator ftIt;
        for (ftIt = ftGroupIt->second.begin(); ftIt != ftGroupIt->second.end(); ++ftIt) {
          featureToSpectrumCluster[*ftIt].push_back(clusterIdx);
        }
      }
    }
  }
}

void FeatureGroups::printFeatureGroups(
    const std::string& featureGroupsOutFile,
    std::vector<DinosaurFeatureList>& allFeatures,
    const std::map<FeatureId, std::vector<int> >& featureToSpectrumCluster,
    std::map<int, std::vector<DinosaurFeature> >& spectrumClusterToConsensusFeatures) {  
  std::ofstream dataStream(featureGroupsOutFile.c_str(), ios::out);
  dataStream.precision(9);
  
  std::vector<std::vector<FeatureIdMatch> >::iterator ftGroupIt;
  size_t featureGroupIdx = 1u;
  for (ftGroupIt = featureGroups_.begin(); ftGroupIt != featureGroups_.end(); ++ftGroupIt) {
    if (ftGroupIt->empty()) continue;
    
    SimilarityMatrix<FeatureId> simMatrix;
    std::set<FeatureId> featureIds;
    
    std::vector<FeatureIdMatch>::const_iterator ftMatchIt;
    for (ftMatchIt = ftGroupIt->begin(); ftMatchIt != ftGroupIt->end(); ++ftMatchIt) {
      featureIds.insert(ftMatchIt->queryFeatureId);
      featureIds.insert(ftMatchIt->targetFeatureId);
      
      simMatrix.updateLink(ftMatchIt->queryFeatureId, 
                           ftMatchIt->targetFeatureId, 
                           1.0 - ftMatchIt->posteriorErrorProb);
    }
    
    std::map<FeatureId, std::map<int, float> > spectrumClusterLinkPEPs;
    propagateSpectrumClusterIds(featureIds, featureToSpectrumCluster, 
        simMatrix, spectrumClusterLinkPEPs);
          
    float sumPrecMz = 0.0f;
    std::set<int> uniqueFileIdxs;
    std::vector<DinosaurFeature> features;
    std::set<FeatureId>::const_iterator ftIdIt;
    for (ftIdIt = featureIds.begin(); ftIdIt != featureIds.end(); ++ftIdIt) {
      int fileIdx = ftIdIt->fileIdx;
      int featureIdx = ftIdIt->scannr;
      
      const DinosaurFeature& ft = allFeatures.at(fileIdx).at(featureIdx);
      
      if (ft.intensity > 0) { /* do not print placeholders */
        features.push_back(ft);
        uniqueFileIdxs.insert(ft.fileIdx);          
        sumPrecMz += ft.precMz;
      }
    }
    
    size_t numFilesMissing = allFeatures.size() - uniqueFileIdxs.size();
    if (numFilesMissing <= maxMissingValues_ && !features.empty()) {        
      DinosaurFeature consensusFeature = features.at(0);
      consensusFeature.precMz = sumPrecMz / features.size();
      consensusFeature.featureIdx = featureGroupIdx++;
      consensusFeature.fileIdx = -1;
      
      std::map<int, size_t> spectrumClusterIdxOffsets;
      addConsensusFeatureToSpectrumClusters(
          consensusFeature, spectrumClusterLinkPEPs,
          spectrumClusterToConsensusFeatures, spectrumClusterIdxOffsets);
      
      printFeatureGroup(features, spectrumClusterLinkPEPs, 
          spectrumClusterIdxOffsets, dataStream);
      
      if (Globals::VERB > 2 && featureGroupIdx % 10000 == 0) {
        std::cerr << "Writing feature group " << featureGroupIdx << std::endl;
      }
    }
  }
}

void FeatureGroups::addConsensusFeatureToSpectrumClusters(
    const DinosaurFeature& consensusFeature,
    std::map<FeatureId, std::map<int, float> >& spectrumClusterLinkPEPs,
    std::map<int, std::vector<DinosaurFeature> >& spectrumClusterToConsensusFeatures,
    std::map<int, size_t>& spectrumClusterIdxOffsets) {
  std::map<int, float>::const_iterator specClustIt;
  FeatureId ftId = spectrumClusterLinkPEPs.begin()->first;
  for (specClustIt = spectrumClusterLinkPEPs[ftId].begin(); 
       specClustIt != spectrumClusterLinkPEPs[ftId].end(); 
       ++specClustIt) {
    int specClustIdx = specClustIt->first;
    /* MaRaCluster creates scan numbers by scannr*100+i, this breaks if i>99 */
    if (specClustIdx != 0 &&        
        spectrumClusterToConsensusFeatures[specClustIdx].size() < 99) {
      spectrumClusterToConsensusFeatures[specClustIdx].push_back(consensusFeature);
      spectrumClusterIdxOffsets[specClustIdx] = spectrumClusterToConsensusFeatures[specClustIdx].size();
    }
  }   
}

void FeatureGroups::printFeatureGroup(
    const std::vector<DinosaurFeature>& features,
    std::map<FeatureId, std::map<int, float> >& spectrumClusterLinkPEPs,
    std::map<int, size_t>& spectrumClusterIdxOffsets,
    std::ostream& dataStream) {
  std::map<int, float>::const_iterator specClustIt;
  std::vector<DinosaurFeature>::const_iterator ftIt;        
  for (ftIt = features.begin(); ftIt != features.end(); ++ftIt) {          
    dataStream << ftIt->fileIdx << "\t" 
               << ftIt->precMz << "\t" 
               << ftIt->charge << "\t" 
               << ftIt->rTime << "\t" 
               << ftIt->intensity << "\t";
    
    FeatureId ftId(ftIt->fileIdx, ftIt->featureIdx);
    for (specClustIt = spectrumClusterLinkPEPs[ftId].begin(); 
         specClustIt != spectrumClusterLinkPEPs[ftId].end(); ) {
      int specClustIdx = specClustIt->first;
      dataStream << specClustIdx*100 + spectrumClusterIdxOffsets[specClustIdx] 
                 << ";" << specClustIt->second;
      ++specClustIt;
      if (specClustIt != spectrumClusterLinkPEPs[ftId].end()) {
        dataStream << ",";
      }
    }
    dataStream << std::endl;
  }
  dataStream << std::endl; 
}

void FeatureGroups::propagateSpectrumClusterIds(
    const std::set<FeatureId>& featureIds, 
    const std::map<FeatureId, std::vector<int> >& featureToSpectrumCluster,
    SimilarityMatrix<FeatureId>& simMatrix,
    std::map<FeatureId, std::map<int, float> >& spectrumClusterLinkPEPs) {
  bool foundClusterIdx = false;
  std::vector<int> clusterIdxs;
  for (std::set<FeatureId>::const_iterator ftIt = featureIds.begin(); 
       ftIt != featureIds.end(); ++ftIt) {
    if (featureToSpectrumCluster.find(*ftIt) != featureToSpectrumCluster.end() &&
        !featureToSpectrumCluster.find(*ftIt)->second.empty()) {
      std::vector<int> clusterIdxs(featureToSpectrumCluster.find(*ftIt)->second);
      foundClusterIdx = true;
      std::map<FeatureId, float> similarities;
      simMatrix.computeShortestPathsFromSource(*ftIt, similarities);
      updateLinkPEPs(*ftIt, clusterIdxs, similarities, spectrumClusterLinkPEPs);
    }
  }
  
  if (!foundClusterIdx) {
    /* if no spectra are associated with the feature group, pick one as root
         TODO: replace this by the "most central" node */
    FeatureId rootFeatureId = *(featureIds.begin()); 
    std::vector<int> clusterIdxs(1, 0);
    std::map<FeatureId, float> similarities;
    simMatrix.computeShortestPathsFromSource(rootFeatureId, similarities);
    updateLinkPEPs(rootFeatureId, clusterIdxs, similarities, spectrumClusterLinkPEPs);
  }
}

void FeatureGroups::updateLinkPEPs(
    const FeatureId rootFeatureId,
    const std::vector<int>& clusterIdxs,
    const std::map<FeatureId, float>& similarities,
    std::map<FeatureId, std::map<int, float> >& spectrumClusterLinkPEPs) {
  std::vector<int>::const_iterator clusterIdxIt;
  for (clusterIdxIt = clusterIdxs.begin(); 
       clusterIdxIt != clusterIdxs.end(); ++clusterIdxIt) {
    std::map<FeatureId, float>::const_iterator ftMatchIt;
    for (ftMatchIt = similarities.begin(); ftMatchIt != similarities.end(); ++ftMatchIt) {
      if (rootFeatureId != ftMatchIt->first && 
          (spectrumClusterLinkPEPs[ftMatchIt->first].find(*clusterIdxIt) == spectrumClusterLinkPEPs[ftMatchIt->first].end() || 
           1.0 - ftMatchIt->second < spectrumClusterLinkPEPs[ftMatchIt->first][*clusterIdxIt])) {
        spectrumClusterLinkPEPs[ftMatchIt->first][*clusterIdxIt] = 1.0 - ftMatchIt->second;
      }
    }
    spectrumClusterLinkPEPs[rootFeatureId][*clusterIdxIt] = 0.0;
  }
}

} /* namespace quandenser */
