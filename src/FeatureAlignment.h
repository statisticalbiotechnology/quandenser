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

#ifndef QUANDENSER_FEATUREALIGNMENT_H_
#define QUANDENSER_FEATUREALIGNMENT_H_

#include <vector>
#include <map>
#include <set>

#include <boost/filesystem.hpp>

#include "maracluster/src/SpectrumFileList.h"

#include "AlignRetention.h"
#include "SplineRegression.h"
#include "DinosaurFeatureList.h"
#include "PercolatorAdapter.h"
#include "DinosaurIO.h"

namespace quandenser {

struct FeatureFeatureMatch {
  FeatureFeatureMatch(double s, const DinosaurFeature& qft,
                      const DinosaurFeature& tft, int l) :
    score(s), queryFeature(qft), targetFeature(tft), label(l) {}

  double score;
  DinosaurFeature queryFeature, targetFeature;
  int label;
};

struct FeatureIdxMatch {
  FeatureIdxMatch() : targetFeatureIdx(-1), posteriorErrorProb(1.0),
                      intensity(0.0) {}
  FeatureIdxMatch(int tftIdx, float pep, float i) :
    targetFeatureIdx(tftIdx), posteriorErrorProb(pep), intensity(i) {}

  int targetFeatureIdx;
  float posteriorErrorProb, intensity;
};


bool operator<(const FeatureFeatureMatch& l, const FeatureFeatureMatch& r);

class FeatureAlignment {
 public:
  FeatureAlignment(const std::string& percolatorOutputFileBaseFN,
      std::vector<std::string>& percolatorArgs,
      float ppmTol, float rTimeStdevTol, float decoyOffset,
      float linkPEPThreshold, float linkPEPMbrSearchThreshold,
      int maxFeatureCandidates) :
    percolatorOutputFileBaseFN_(percolatorOutputFileBaseFN),
    percolatorArgs_(percolatorArgs),
    ppmTol_(ppmTol), rTimeStdevTol_(rTimeStdevTol), decoyOffset_(decoyOffset),
    linkPEPThreshold_(linkPEPThreshold),
    linkPEPMbrSearchThreshold_(linkPEPMbrSearchThreshold),
    maxFeatureCandidates_(maxFeatureCandidates) {}

  void matchFeatures(
    const std::vector<std::pair<int, FilePair> >& featureAlignmentQueue,
    const maracluster::SpectrumFileList& fileList,
    AlignRetention& alignRetention,
    std::vector<DinosaurFeatureList>& allFeatures,
    const std::string& addedFeaturesFile, const std::string& tmpFilePrefix);

  std::map<FilePair, std::map<int, FeatureIdxMatch> >& getFeatureMatches() {
    return featureMatches_;
  }

  static size_t loadFromFile(const std::string& ftMatchFile, std::map<int, FeatureIdxMatch>& ftMatchMap) {
    std::vector<std::pair<int, FeatureIdxMatch> > addedFtMatches;
    maracluster::BinaryInterface::read(ftMatchFile, addedFtMatches);
    for (std::pair<int, FeatureIdxMatch>& df : addedFtMatches) {
      ftMatchMap[df.first] = df.second;
    }
    return ftMatchMap.size();
  }

  static void saveToFile(const std::string& ftMatchFile, std::map<int, FeatureIdxMatch>& ftMatchMap, bool append) {
    std::vector<std::pair<int, FeatureIdxMatch> > featureMatchPairs(
        ftMatchMap.begin(), ftMatchMap.end());
    maracluster::BinaryInterface::write(featureMatchPairs, ftMatchFile, append);
  }

 protected:
  std::string percolatorOutputFileBaseFN_;
  std::vector<std::string>& percolatorArgs_;

  float ppmTol_, rTimeStdevTol_;
  float decoyOffset_;
  float linkPEPThreshold_, linkPEPMbrSearchThreshold_;
  int maxFeatureCandidates_;

  std::map<FilePair, std::map<int, FeatureIdxMatch> > featureMatches_;

  static const std::vector<std::pair<std::string, double> > kLinkFeatureNames;

  void getFeatureMap(FilePair& filePair,
    const std::string& targetMzMLFile,
    SplineRegression& alignment,
    DinosaurFeatureList& featuresQueryRun,
    DinosaurFeatureList& featuresTargetRun,
		const std::string& addedFeaturesFile,
    const std::string& tmpFilePrefix);

  void insertPlaceholderFeatures(
    const DinosaurFeatureList& featuresQueryRun,
    const std::vector<double>& predictedRTimesTargetRun,
    DinosaurFeatureList& featuresTargetRun,
    std::map<int, FeatureIdxMatch>& precursorLinks,
    const std::string& addedFeaturesFile);

  void matchFeatures(
    const std::string& percolatorOutputFile,
    const DinosaurFeatureList& featuresQueryRun,
    const DinosaurFeatureList& featuresTargetRun,
    const std::vector<double>& predictedRTimesTargetRun,
    float rTimeStdev);

  void addFeatureLinks(const std::string& percolatorOutputFile,
    DinosaurFeatureList& candidateFeaturesTargetRun,
    DinosaurFeatureList& featuresTargetRun,
    std::map<int, FeatureIdxMatch>& precursorLinks,
    const std::string& addedFeaturesFile);

  void findMatches(int label, double rTimeTol,
    const DinosaurFeatureList& featuresTargetRun,
    const DinosaurFeature& queryFeature,
    std::vector<FeatureFeatureMatch>& topTargets);

  inline float getPpmDiff(float precMz1, float precMz2) {
    return std::abs(precMz1 - precMz2) / precMz1 * 1e6;
  }
  double getXicRateDiff(const DinosaurFeature& queryFeature,
    const DinosaurFeature& targetFeature);
  double getXicRate(const double intensity,
    const double rtStart, const double rtEnd);

  /* mbr: matches-between-runs */
  void mbrMatchFeatures(
    const std::map<int, FeatureIdxMatch>& precursorLinks,
    const std::string& percolatorOutputFile,
    const std::string& targetMzMLFile,
    const DinosaurFeatureList& featuresQueryRun,
    const DinosaurFeatureList& featuresTargetRun,
    const std::vector<double>& predictedRTimesTargetRun,
    float rTimeStdev,
    DinosaurFeatureList& candidateFeaturesTargetRun);

  void processDinosaurTargets(
    const std::map<int, FeatureIdxMatch>& precursorLinks,
    const std::string& dinosaurTargetOutputFile,
    const std::string& percolatorOutputFile,
    const DinosaurFeatureList& featuresQueryRun,
    const DinosaurFeatureList& featuresTargetRun,
    const std::vector<double>& predictedRTimesTargetRun,
    DinosaurFeatureList& candidateFeaturesTargetRun);

  void addLinkPsm(int label,
    const DinosaurFeature& queryFeature,
    const DinosaurFeature& targetFeature,
    PercolatorAdapter& percolatorAdapter);

  std::string convertFeatureToPsmId(const DinosaurFeature& feature);
  DinosaurFeature parsePsmIdAsFeature(const std::string& specId);

  /* start negative candidate feature indices at -2 */
  static inline int getCandidateIdx(int size) { return -1*(size + 2); }
  static inline size_t getCandidatePos(int featureIdx) {
    return static_cast<size_t>(featureIdx * -1 - 2);
  }

};

} /* namespace quandenser */

#endif /* QUANDENSER_FEATUREALIGNMENT_H_ */
