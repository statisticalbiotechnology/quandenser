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

#ifndef QUANDENSER_MARACLUSTERIO_H_
#define QUANDENSER_MARACLUSTERIO_H_

#include <vector>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "percolator/src/DataSet.h" /* contains TabReader class */

#include "SpectrumFiles.h"
#include "DinosaurFeature.h"
#include "DinosaurFeatureList.h"
#include "SpectrumToPrecursorMap.h"
#include "AlignRetention.h"

namespace quandenser {

class MaRaClusterIO {
 public:
  static void parseClustersForRTimePairs(std::istream& dataStream, 
    maracluster::SpectrumFileList& fileList, 
    SpectrumToPrecursorMap& spectrumToPrecursorMap, 
    RTimePairs& rTimePairs);
  static void parseClustersForConsensusMerge(std::istream& dataStream, 
    maracluster::SpectrumFileList& fileList, 
    SpectrumToPrecursorMap& spectrumToPrecursorMap,
    std::map<maracluster::ScanId, std::vector<int> >& featureToSpectrumClusterIdxMap);
  
  static void setPrecursorTolerance(int ppm) {
    precursorTolerancePpm_ = ppm;
  }
  
 protected:  
  static void generateRTimePairs(DinosaurFeatureList& dinosaurFeatures, RTimePairs& rTimePairs);
  static maracluster::ScanId parseClusterRow(std::string& line, maracluster::SpectrumFileList& fileList);
  static float precursorTolerancePpm_;
  
};

} /* namespace quandenser */

#endif /* QUANDENSER_MARACLUSTERADAPTER_H_ */
