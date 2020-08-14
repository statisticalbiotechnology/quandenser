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

#ifndef QUANDENSER_QUANDENSER_H_
#define QUANDENSER_QUANDENSER_H_

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include <boost/filesystem.hpp>
#include <boost/asio.hpp>
#include <boost/functional/hash_fwd.hpp>

#include "maracluster/src/Globals.h"
#include "maracluster/src/SpectrumFileList.h"

#include "Globals.h"
#include "Option.h"
#include "Version.h"

#include "MaRaClusterAdapter.h"
#include "SpectrumFiles.h"
#include "DinosaurIO.h"
#include "AlignRetention.h"
#include "SplineRegression.h"
#include "FeatureAlignment.h"
#include "FeatureGroups.h"
#include "MaRaClusterIO.h"

namespace quandenser {


class Quandenser {
 public:
  Quandenser();
  ~Quandenser() {};

  bool parseOptions(int argc, char **argv);

  int run();

 protected:
  std::string greeter();
  std::string extendedGreeter(time_t& startTime);

  std::string call_;
  std::string spectrumBatchFileFN_;

  std::string outputFolder_;
  std::string outputSpectrumFile_;
  std::string fnPrefix_;
  int seed_;
  unsigned int numThreads_;
  int maxMissingValues_;
  float intensityScoreThreshold_;

  float maraclusterPpmTol_, alignPpmTol_, alignRTimeStdevTol_;
  float decoyOffset_;
  float linkPEPThreshold_, linkPEPMbrSearchThreshold_;
  int maxFeatureCandidates_;

  int partial1Dinosaur_;
  std::string partial2MaRaCluster_;
  int partial3MatchRound_;
  std::string partial4Consensus_;
	bool useTempFiles_;

  int runMaRaCluster(const std::string& maRaClusterSubFolder,
    const maracluster::SpectrumFileList& fileList,
    std::vector<DinosaurFeatureList>& allFeatures,
    std::string& clusterFilePath,
    const std::string& spectrumToPrecursorFile,
    const std::string& tmpFilePrefixAlign,
    const std::string& mode);

  void loadAllFeatures(const std::string& tmpFilePrefixAlign,
    std::vector<DinosaurFeatureList>& allFeatures);
  void unloadAllFeatures(std::vector<DinosaurFeatureList>& allFeatures);
  
  int saveAlignmentQueue(
    std::vector<std::pair<int, FilePair> >& featureAlignmentQueue,
    const std::string& featureAlignQueueFile,
    maracluster::SpectrumFileList& fileList);
  
  void updateMatchesTmpFile(
    const boost::filesystem::path& matchesTmpFileBaserefixAlign, 
    const FilePair& filePair, 
    const DinosaurFeatureList& currentFeatures, 
    const DinosaurFeatureList& addedFeatures,
    const size_t originalNumFeaturesFromFile);
    
  std::string getFeatureFN(const boost::filesystem::path& basePath,
                           size_t fileIdx);
  std::string getAddedFeaturesFN(const boost::filesystem::path& basePath,
                           size_t fileIdx1, size_t fileIdx2);
  
  int createDirectory(boost::filesystem::path dirPath);
  
  // google analytics
  static bool parseUrl(std::string url, std::string* host, std::string* path);
  static void httpRequest(const std::string& url, const std::string& data);
  static void postToAnalytics(const std::string& appName);

  std::vector<std::string> maraclusterArgs_;
  std::vector<std::string> percolatorArgs_;
};

} /* namespace quandenser */

#endif /* QUANDENSER_QUANDENSER_H_ */
