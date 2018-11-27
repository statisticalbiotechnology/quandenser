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

#include "MaRaClusterIO.h"

namespace quandenser {

float MaRaClusterIO::precursorTolerancePpm_ = 20.0f;

void MaRaClusterIO::parseClustersForRTimePairs(std::istream& dataStream, 
    maracluster::SpectrumFileList& fileList, 
    SpectrumToPrecursorMap& spectrumToPrecursorMap, 
    RTimePairs& rTimePairs) {
  std::string clusterLine;
  
  DinosaurFeatureList dinosaurFeatures;
  size_t lineNr = 1;
  size_t clusterIdx = 1;
  while (std::getline(dataStream, clusterLine)) {
    if (lineNr % 100000 == 0 && Globals::VERB > 1) {
      std::cerr << "Processing line " << lineNr << std::endl;
    }
    if (clusterLine.size() > 0) {
      maracluster::ScanId spectrumId = parseClusterRow(clusterLine, fileList);
      DinosaurFeatureList::const_iterator ftIt;
      for (ftIt = spectrumToPrecursorMap.getBegin(spectrumId); 
           ftIt != spectrumToPrecursorMap.getEnd(spectrumId); ++ftIt) {
        dinosaurFeatures.push_back(*ftIt);
      }
    } else {
      generateRTimePairs(dinosaurFeatures, rTimePairs);
      dinosaurFeatures.clear();
      ++clusterIdx;
    }
    
    ++lineNr;
  }
  
  generateRTimePairs(dinosaurFeatures, rTimePairs);
  
  if (Globals::VERB > 2) {
    for (int i = 0; i < fileList.size(); ++i) {
      for (int j = i+1; j < fileList.size(); ++j) {
        std::cerr << "File pair (" << i << "," << j << "): " << 
            rTimePairs[FilePair(i, j)].size() << " retention time pairs." << std::endl;
      }
    }
  }
}

void MaRaClusterIO::generateRTimePairs(DinosaurFeatureList& dinosaurFeatures, 
    RTimePairs& rTimePairs) {
  for (DinosaurFeatureList::const_iterator it = dinosaurFeatures.begin(); 
       it != dinosaurFeatures.end(); ++it) {
    float lowerPrecMzBound = it->precMz*(1.0 - precursorTolerancePpm_ * 1e-6); 
    float upperPrecMzBound = it->precMz*(1.0 + precursorTolerancePpm_ * 1e-6);
    for (DinosaurFeatureList::const_iterator it2 = it + 1; 
         it2 != dinosaurFeatures.end(); ++it2) {
      if (it->charge == it2->charge && lowerPrecMzBound <= it2->precMz
            && upperPrecMzBound >= it2->precMz) {
        int fileIdx1 = it->fileIdx;
        float rTime1 = it->rTime;
        int fileIdx2 = it2->fileIdx;
        float rTime2 = it2->rTime;
        
        if (fileIdx1 == fileIdx2) {
          continue;
        } else if (fileIdx1 > fileIdx2) {
          std::swap(fileIdx1, fileIdx2);
          std::swap(rTime1, rTime2);
        }
        rTimePairs[FilePair(fileIdx1, fileIdx2)].push_back(RTimePair(rTime1, rTime2));
      }
    }
  }
}

maracluster::ScanId MaRaClusterIO::parseClusterRow(std::string& line, 
    maracluster::SpectrumFileList& fileList) {
  TabReader reader(line);
  
  std::string filepath = reader.readString();
  int scannr = reader.readInt();
  
  return fileList.getScanId(filepath, scannr);
}

void MaRaClusterIO::parseClustersForConsensusMerge(std::istream& dataStream, 
    maracluster::SpectrumFileList& fileList, 
    SpectrumToPrecursorMap& spectrumToPrecursorMap,
    std::map<maracluster::ScanId, std::vector<int> >& featureToSpectrumCluster) {
  std::string clusterLine;
  size_t clusterIdx = 1;
  if (Globals::VERB > 1) {
    std::cerr << "Mapping features to spectra." << std::endl;
  }
  while (std::getline(dataStream, clusterLine)) {
    if (clusterLine.size() > 0) {
      maracluster::ScanId spectrumId = parseClusterRow(clusterLine, fileList);
      DinosaurFeatureList::const_iterator ftIt;
      for (ftIt = spectrumToPrecursorMap.getBegin(spectrumId); 
           ftIt != spectrumToPrecursorMap.getEnd(spectrumId); ++ftIt) {
        featureToSpectrumCluster[maracluster::ScanId(ftIt->fileIdx, ftIt->featureIdx)].push_back(clusterIdx);
      }
    } else {
      ++clusterIdx;
    }
  }
}

} /* namespace quandenser */
