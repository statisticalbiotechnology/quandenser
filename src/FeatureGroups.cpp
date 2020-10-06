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

#include <unistd.h>
#include <sys/resource.h>

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( )
{
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);
}

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
  std::vector<std::set<size_t> > linksToObservedFeature(numFiles);
  std::map<size_t, size_t> groupIdMap;
  
  if (!tmpFilePrefixGroup.empty()) {
    for (size_t i = 0; i < numFiles; ++i) {
      std::string fileName = tmpFilePrefixGroup + boost::lexical_cast<std::string>(i) + ".dat";
      remove(fileName.c_str());
    }
  }
  
  std::cerr << "Mem usage: " << getCurrentRSS() << std::endl;
  
  std::vector<std::pair<int, FilePair> >::const_iterator filePairIt;
  int filePairCnt = 0;
  for (filePairIt = featureAlignmentQueue.begin(); 
       filePairIt != featureAlignmentQueue.begin() + (featureAlignmentQueue.size() / 2); 
       ++filePairIt, ++filePairCnt) {
    FilePair filePair = filePairIt->second;
    int queryFileIdx = filePair.fileIdx1;
    int targetFileIdx = filePair.fileIdx2;
    
    FeatureToGroupMap& queryFeatureIdToGroupId = featureIdToGroupId.at(queryFileIdx);
    FeatureToGroupMap& targetFeatureIdToGroupId = featureIdToGroupId.at(targetFileIdx);
    
    if (!tmpFilePrefixGroup.empty()) {
      queryFeatureIdToGroupId.loadFromFile(tmpFilePrefixGroup + boost::lexical_cast<std::string>(queryFileIdx) + ".dat");
      targetFeatureIdToGroupId.loadFromFile(tmpFilePrefixGroup + boost::lexical_cast<std::string>(targetFileIdx) + ".dat");
    }
    
    if (Globals::VERB > 2) {
      std::cerr << "Processing pair " << filePairCnt+1 << "/" 
          << featureAlignmentQueue.size() / 2 << std::endl;
    }
    
    bool isRevFilePair = false;
    processFeatureAlignment(filePair, featureMatches[filePair], 
        tmpFilePrefixAlign, queryFeatureIdToGroupId, targetFeatureIdToGroupId,
        linksToObservedFeature.at(queryFileIdx), linksToObservedFeature.at(targetFileIdx), 
        groupIdMap, isRevFilePair);
    
    FilePair revFilePair = filePair.getRevFilePair();
    isRevFilePair = true;
    processFeatureAlignment(revFilePair, featureMatches[revFilePair], 
        tmpFilePrefixAlign, queryFeatureIdToGroupId, targetFeatureIdToGroupId,
        linksToObservedFeature.at(queryFileIdx), linksToObservedFeature.at(targetFileIdx), 
        groupIdMap, isRevFilePair);
    
    if (!tmpFilePrefixGroup.empty()) {
      bool append = false;
      queryFeatureIdToGroupId.saveToFile(tmpFilePrefixGroup + boost::lexical_cast<std::string>(queryFileIdx) + ".dat", append);
      targetFeatureIdToGroupId.saveToFile(tmpFilePrefixGroup + boost::lexical_cast<std::string>(targetFileIdx) + ".dat", append);
      
      queryFeatureIdToGroupId.clear();
      targetFeatureIdToGroupId.clear();
    }
    
    linksToObservedFeature.at(queryFileIdx).clear();
    
    std::cerr << "Num clusters: " << featureGroups_.size() << std::endl;
    std::cerr << "Num clusters allocated: " << featureGroups_.capacity() << std::endl;
    std::vector<std::vector<FeatureIdMatch> >::iterator ftGroupIt;
    size_t numLinks = 0;
    size_t numLinksAllocated = 0;
    for (ftGroupIt = featureGroups_.begin(); ftGroupIt != featureGroups_.end(); ++ftGroupIt) {
      numLinks += ftGroupIt->size();
      numLinksAllocated += ftGroupIt->capacity();
    }
    std::cerr << "Num links: " << numLinks << std::endl;
    std::cerr << "Num links allocated: " << numLinksAllocated << std::endl;
    std::cerr << "Num cluster id mappings: " << groupIdMap.size() << std::endl;
    std::cerr << "Mem usage: " << getCurrentRSS() << std::endl;
  }
  
  if (!tmpFilePrefixGroup.empty()) {
    for (size_t i = 0; i < numFiles; ++i) {
      std::string fileName = tmpFilePrefixGroup + boost::lexical_cast<std::string>(i) + ".dat";
      remove(fileName.c_str());
    }
  }
}

void FeatureGroups::processFeatureAlignment(FilePair& filePair, 
    std::map<int, FeatureIdxMatch>& filePairFeatureMatches,
    const std::string& tmpFilePrefixAlign,
    FeatureToGroupMap& queryFeatureIdToGroupId,
    FeatureToGroupMap& targetFeatureIdToGroupId,
    std::set<size_t>& queryLinksToObservedFeature,
    std::set<size_t>& targetLinksToObservedFeature,
    std::map<size_t, size_t>& groupIdMap,
    bool isRevFilePair) {
  int queryFileIdx = filePair.fileIdx1;
  int targetFileIdx = filePair.fileIdx2;
  
  if (!tmpFilePrefixAlign.empty()) {
    std::string matchesFileName = FeatureAlignment::getAddedFeaturesFN(
      tmpFilePrefixAlign + "/matches", queryFileIdx, targetFileIdx);
    FeatureAlignment::loadFromFile(matchesFileName, filePairFeatureMatches);
  }
  
  if (Globals::VERB > 2) {
    std::cerr << "  " << queryFileIdx << "->" << targetFileIdx << " (" 
        << filePairFeatureMatches.size() << " links)" << std::endl;
  }
  
  if (isRevFilePair) std::swap(queryFileIdx, targetFileIdx);
  
  std::map<int, FeatureIdxMatch>::const_iterator ftIdxMatchIt;
  for (ftIdxMatchIt = filePairFeatureMatches.begin(); 
       ftIdxMatchIt != filePairFeatureMatches.end();
       ++ftIdxMatchIt) {
    int queryFeatureIdx = ftIdxMatchIt->first;
    int targetFeatureIdx = ftIdxMatchIt->second.targetFeatureIdx;
    float posteriorErrorProb = ftIdxMatchIt->second.posteriorErrorProb;
    
    if (isRevFilePair) std::swap(queryFeatureIdx, targetFeatureIdx);
    
    FeatureId queryFeatureId(queryFileIdx, queryFeatureIdx);
    FeatureId targetFeatureId(targetFileIdx, targetFeatureIdx);
    
    /* placeholder links have a PEP of 1.0, we do not have to follow these
       if no observed feature is present further down the alignment tree */    
    if (posteriorErrorProb == 1.0) {
      if (isRevFilePair && queryLinksToObservedFeature.find(queryFeatureIdx) == queryLinksToObservedFeature.end()) {
        continue;
      }
      /* PEP=1.0 would destroy the linkPEP score, we deal with this downstream */
      posteriorErrorProb = 0.0;
    }
    targetLinksToObservedFeature.insert(targetFeatureIdx);
    
    addToFeatureGroup(queryFeatureId, targetFeatureId, posteriorErrorProb,
        queryFeatureIdToGroupId, targetFeatureIdToGroupId, groupIdMap);
  }
  
  filePairFeatureMatches.clear();
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
  std::vector<FeatureIdMatch>().swap(featureGroups_[mergeInClusterIdx]);
}

void FeatureGroups::filterConsensusFeatures(
    const std::vector<DinosaurFeatureList>& allFeatures,
    std::map<FeatureId, std::vector<int> >& featureToSpectrumCluster) {
  if (Globals::VERB > 2) {
    std::cerr << "Applying intensity score filter for consensus spectrum "
              << "to feature group assignments." << std::endl;
    std::cerr << "  Mem usage: " << getCurrentRSS() << std::endl;
  }
  
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
        int vectorIdx = allFeatures.at(fileIdx).getVectorIdx(featureIdx);
        if (vectorIdx >= 0) {
          const DinosaurFeature& ft = allFeatures.at(fileIdx).at(vectorIdx);
          
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
  
  if (Globals::VERB > 2) {
    std::cerr << "Finished applying intensity score filter." << std::endl;
    std::cerr << "  Mem usage: " << getCurrentRSS() << std::endl;
  }
}

void FeatureGroups::printFeatureGroups(
    const std::string& featureGroupsOutFile,
    const std::vector<DinosaurFeatureList>& allFeatures,
    const std::map<FeatureId, std::vector<int> >& featureToSpectrumCluster,
    std::map<int, std::vector<DinosaurFeature> >& spectrumClusterToConsensusFeatures) {  
  if (Globals::VERB > 2) {
    std::cerr << "Writing feature groups to output file." << std::endl;
    std::cerr << "  Mem usage: " << getCurrentRSS() << std::endl;
  }
  
  std::ofstream dataStream(featureGroupsOutFile.c_str(), ios::out);
  dataStream.precision(9);
  
  size_t featureGroupIdx = 1u;
#pragma omp parallel for schedule(dynamic, 1000) 
  for (size_t processedFeatureGroups = 0u; processedFeatureGroups < featureGroups_.size(); ++processedFeatureGroups) {
    std::vector<FeatureIdMatch>& ftGroup = featureGroups_.at(processedFeatureGroups);
    if (Globals::VERB > 2 && processedFeatureGroups % 10000 == 0) {
      std::cerr << "  Processing feature group " << processedFeatureGroups
                << " / " << featureGroups_.size() << std::endl;
    }
    
    if (ftGroup.empty()) continue;
    
    SimilarityMatrix<FeatureId> simMatrix;
    std::set<FeatureId> featureIds;
    
    std::vector<FeatureIdMatch>::const_iterator ftMatchIt;
    for (ftMatchIt = ftGroup.begin(); ftMatchIt != ftGroup.end(); ++ftMatchIt) {
      featureIds.insert(ftMatchIt->queryFeatureId);
      featureIds.insert(ftMatchIt->targetFeatureId);
      
      simMatrix.updateLink(ftMatchIt->queryFeatureId, 
                           ftMatchIt->targetFeatureId, 
                           1.0 - ftMatchIt->posteriorErrorProb);
    }
          
    float sumPrecMz = 0.0f;
    std::set<int> uniqueFileIdxs;
    std::vector<DinosaurFeature> features;
    std::set<FeatureId>::const_iterator ftIdIt;
    for (ftIdIt = featureIds.begin(); ftIdIt != featureIds.end(); ++ftIdIt) {
      int fileIdx = ftIdIt->fileIdx;
      int featureIdx = ftIdIt->scannr;
      
      int vectorIdx = allFeatures.at(fileIdx).getVectorIdx(featureIdx);
      if (vectorIdx >= 0) { /* do not print placeholders */
        const DinosaurFeature& ft = allFeatures.at(fileIdx).at(vectorIdx);
        features.push_back(ft);
        uniqueFileIdxs.insert(ft.fileIdx);          
        sumPrecMz += ft.precMz;
      }
    }
    
    size_t numFilesMissing = allFeatures.size() - uniqueFileIdxs.size();
    if (numFilesMissing <= maxMissingValues_ && !features.empty()) {        
      DinosaurFeature consensusFeature = features.at(0);
      consensusFeature.precMz = sumPrecMz / features.size();
    #pragma omp critical(feature_idx)
      {
        consensusFeature.featureIdx = featureGroupIdx++;
        if (Globals::VERB > 2 && featureGroupIdx % 10000 == 0) {
          std::cerr << "  Writing feature group " << featureGroupIdx << std::endl;
        }
      }
      consensusFeature.fileIdx = -1;
      
      std::map<int, std::map<FeatureId, float> > spectrumClusterLinkPEPs;
      propagateSpectrumClusterIds(featureIds, featureToSpectrumCluster, 
          simMatrix, spectrumClusterLinkPEPs);
      
      std::map<int, size_t> spectrumClusterIdxOffsets;
    #pragma omp critical(assign_features)
      {
        addConsensusFeatureToSpectrumClusters(
            consensusFeature, spectrumClusterLinkPEPs,
            spectrumClusterToConsensusFeatures, spectrumClusterIdxOffsets);
      }
    
    #pragma omp critical(print_features)
      {
        printFeatureGroup(features, spectrumClusterLinkPEPs, 
            spectrumClusterIdxOffsets, dataStream);
      }
    }
  }
  
  if (Globals::VERB > 2) {
    std::cerr << "Finished writing feature groups." << std::endl;
    std::cerr << "  Mem usage: " << getCurrentRSS() << std::endl;
  }
}

void FeatureGroups::addConsensusFeatureToSpectrumClusters(
    const DinosaurFeature& consensusFeature,
    std::map<int, std::map<FeatureId, float> >& spectrumClusterLinkPEPs,
    std::map<int, std::vector<DinosaurFeature> >& spectrumClusterToConsensusFeatures,
    std::map<int, size_t>& spectrumClusterIdxOffsets) {
  std::map<int, std::map<FeatureId, float> >::const_iterator specClustIt;
  for (specClustIt = spectrumClusterLinkPEPs.begin(); 
       specClustIt != spectrumClusterLinkPEPs.end(); 
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
    std::map<int, std::map<FeatureId, float> >& spectrumClusterLinkPEPs,
    std::map<int, size_t>& spectrumClusterIdxOffsets,
    std::ostream& dataStream) {
  std::map<int, std::map<FeatureId, float> >::iterator specClustIt;
  std::vector<DinosaurFeature>::const_iterator ftIt;        
  for (ftIt = features.begin(); ftIt != features.end(); ++ftIt) {          
    dataStream << ftIt->fileIdx << "\t" 
               << ftIt->precMz << "\t" 
               << ftIt->charge << "\t" 
               << ftIt->rTime << "\t" 
               << ftIt->intensity << "\t";
    
    FeatureId ftId(ftIt->fileIdx, ftIt->featureIdx);
    for (specClustIt = spectrumClusterLinkPEPs.begin(); 
         specClustIt != spectrumClusterLinkPEPs.end(); ) {
      int specClustIdx = specClustIt->first;
      float postErrProb = specClustIt->second[ftId];
      if (postErrProb > 0) {
        postErrProb = 1.0 - postErrProb;
      }
      dataStream << specClustIdx*100 + spectrumClusterIdxOffsets[specClustIdx] 
                 << ";" << postErrProb;
      ++specClustIt;
      if (specClustIt != spectrumClusterLinkPEPs.end()) {
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
    std::map<int, std::map<FeatureId, float> >& spectrumClusterLinkPEPs) {
  bool foundClusterIdx = false;
  std::set<int> featureGroupClusterIdxs;
  
  if (featureIds.size() > 200) {
    std::cerr << "Start" << std::endl;
  }
  for (std::set<FeatureId>::const_iterator ftIt = featureIds.begin(); 
       ftIt != featureIds.end(); ++ftIt) {
    if (featureToSpectrumCluster.find(*ftIt) != featureToSpectrumCluster.end() &&
        !featureToSpectrumCluster.find(*ftIt)->second.empty()) {
      std::vector<int>::const_iterator clusterIt;
      const std::vector<int>& featureClusterIdxs = featureToSpectrumCluster.find(*ftIt)->second;
      for (clusterIt = featureClusterIdxs.begin(); clusterIt != featureClusterIdxs.end(); ++clusterIt) {
        foundClusterIdx = true;
        spectrumClusterLinkPEPs[*clusterIt][*ftIt] = 0.0;
      }
    }
  }
  
  if (featureIds.size() > 200) {
    std::cerr << "Finish " << featureIds.size() << " " << spectrumClusterLinkPEPs.size() << std::endl;
  }
  
  if (!foundClusterIdx) {
    /* if no spectra are associated with the feature group, pick one as root
         TODO: replace this by the "most central" node */
    FeatureId rootFeatureId = *(featureIds.begin());
    spectrumClusterLinkPEPs[0][rootFeatureId] = 0.0;
  }
  
  std::map<int, std::map<FeatureId, float> >::iterator clusterIt;
  for (clusterIt = spectrumClusterLinkPEPs.begin(); clusterIt != spectrumClusterLinkPEPs.end(); ++clusterIt) {
    simMatrix.computeShortestPathsFromSource(clusterIt->second);
  }
}

} /* namespace quandenser */
