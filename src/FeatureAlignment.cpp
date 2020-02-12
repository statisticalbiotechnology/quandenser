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

#include "FeatureAlignment.h"

#include <boost/lexical_cast.hpp>
#include <boost/assign/list_of.hpp>

namespace quandenser {

const std::vector<std::pair<std::string, double> > 
FeatureAlignment::kLinkFeatureNames =
  boost::assign::list_of(make_pair("ppmDiff", -3.0))
    (make_pair("rtDiff",-1.0))
    (make_pair("precMz", 0.0))
    (make_pair("rTime", 0.0))
    /* (make_pair("xicRateDiff", -0.5)) */
    (make_pair("queryIsPlaceHolder", 0.0))
    (make_pair("targetIsPlaceHolder", 0.0))
    (make_pair("charge1", 0.0))
    (make_pair("charge2", 0.0))
    (make_pair("charge3", 0.0))
    (make_pair("charge4plus", 0.0));
    
bool operator<(const FeatureFeatureMatch& l, const FeatureFeatureMatch& r) {
  return l.score < r.score;
}
    
void FeatureAlignment::matchFeatures(
    const std::vector<std::pair<int, FilePair> >& featureAlignmentQueue,
    const maracluster::SpectrumFileList& fileList,
    AlignRetention& alignRetention,
    std::vector<DinosaurFeatureList>& allFeatures,
    const std::string& tmpFilePrefix) {  
  std::vector<std::pair<int, FilePair> >::const_iterator filePairIt;
  size_t alignmentCnt = 1u;
  for (filePairIt = featureAlignmentQueue.begin(); 
        filePairIt != featureAlignmentQueue.end(); ++filePairIt, ++alignmentCnt) {
    FilePair filePair = filePairIt->second;
    if (Globals::VERB > 1) {
      std::cerr << "Matching features " << filePair.fileIdx1 
                << "->" << filePair.fileIdx2 
                << " (" << alignmentCnt << "/" 
                << featureAlignmentQueue.size() << ")" << std::endl;
    }
    std::string targetMzMLFile = fileList.getFilePath(filePair.fileIdx2);
    getFeatureMap(filePair, targetMzMLFile, 
        alignRetention.getAlignment(filePair), 
        allFeatures.at(filePair.fileIdx1), 
        allFeatures.at(filePair.fileIdx2),
        tmpFilePrefix);
  }
}

void FeatureAlignment::getFeatureMap(FilePair& filePair,
    const std::string& targetMzMLFile,
    SplineRegression& alignment, 
    DinosaurFeatureList& featuresQueryRun, 
    DinosaurFeatureList& featuresTargetRun,
    const std::string& tmpFilePrefix) {
  if (!tmpFilePrefix.empty()) {
    std::string fileName = tmpFilePrefix + "/features."  + boost::lexical_cast<std::string>(filePair.fileIdx1) + ".dat";
    featuresQueryRun.loadFromFile(fileName);
    
    bool withIdxMap = true;
    fileName = tmpFilePrefix + "/features."  + boost::lexical_cast<std::string>(filePair.fileIdx2) + ".dat";
    featuresTargetRun.loadFromFile(fileName, withIdxMap);
  }
  
  std::vector<double> rTimesQueryRun, predictedRTimesTargetRun;
  DinosaurFeatureList::const_iterator ftIt;
  for (ftIt = featuresQueryRun.begin(); ftIt != featuresQueryRun.end(); ++ftIt) {
    rTimesQueryRun.push_back(ftIt->rTime);
  }
  alignment.predict(rTimesQueryRun, predictedRTimesTargetRun);
  float rTimeStdev = alignment.getRmse();
  DinosaurFeatureList candidateFeaturesTargetRun;
  
  std::string percolatorOutputFile = percolatorOutputFileBaseFN_ + "/link_" + 
    boost::lexical_cast<std::string>(filePair.fileIdx1) + "_to_" + 
    boost::lexical_cast<std::string>(filePair.fileIdx2) + ".psms";
  if (Globals::fileIsEmpty(percolatorOutputFile)) {
    matchFeatures(percolatorOutputFile, featuresQueryRun, featuresTargetRun, 
                 predictedRTimesTargetRun, rTimeStdev);
  }
  addFeatureLinks(percolatorOutputFile, candidateFeaturesTargetRun, 
                  featuresTargetRun, featureMatches_[filePair]);
  
  if (linkPEPMbrSearchThreshold_ < 1.0) {
    std::string percolatorMbrOutputFile = percolatorOutputFileBaseFN_ + 
      "/search_and_link_" + 
      boost::lexical_cast<std::string>(filePair.fileIdx1) + "_to_" + 
      boost::lexical_cast<std::string>(filePair.fileIdx2) + ".psms.pout";
    mbrMatchFeatures(featureMatches_[filePair], 
        percolatorMbrOutputFile, targetMzMLFile, 
        featuresQueryRun, featuresTargetRun, predictedRTimesTargetRun, 
        rTimeStdev, candidateFeaturesTargetRun);
    
    addFeatureLinks(percolatorMbrOutputFile, candidateFeaturesTargetRun, 
        featuresTargetRun, featureMatches_[filePair]);
  }
  
  insertPlaceholderFeatures(featuresQueryRun, predictedRTimesTargetRun, 
      featuresTargetRun, featureMatches_[filePair]);
  
  featuresTargetRun.sortByPrecMz();
  
  if (!tmpFilePrefix.empty()) {
    std::string fileName = tmpFilePrefix + "/features."  + boost::lexical_cast<std::string>(filePair.fileIdx2) + ".dat";
    bool append = false;
    featuresTargetRun.saveToFile(fileName, append);
    
    featuresQueryRun.clear();
    featuresTargetRun.clear();
    
    std::string matchesFileName = tmpFilePrefix + "/matches."  + 
        boost::lexical_cast<std::string>(filePair.fileIdx1) + "." + 
        boost::lexical_cast<std::string>(filePair.fileIdx2) + ".dat";
    saveToFile(matchesFileName, featureMatches_[filePair], append);
    featureMatches_[filePair].clear();
  }
}

void FeatureAlignment::insertPlaceholderFeatures(
    const DinosaurFeatureList& featuresQueryRun, 
    const std::vector<double>& predictedRTimesTargetRun,
    DinosaurFeatureList& featuresTargetRun,
    std::map<int, FeatureIdxMatch>& precursorLinks) {
  DinosaurFeatureList::const_iterator queryIt;
  std::vector<double>::const_iterator predRTimeIt;
  int targetFileIdx = featuresTargetRun.begin()->fileIdx;
  size_t placeholdersAdded = 0u;
  for (queryIt = featuresQueryRun.begin(), predRTimeIt = predictedRTimesTargetRun.begin(); 
       queryIt != featuresQueryRun.end(); ++queryIt, ++predRTimeIt) {
    if (precursorLinks[queryIt->featureIdx].posteriorErrorProb > linkPEPThreshold_) {
      DinosaurFeature placeholderFt = *queryIt;
      placeholderFt.fileIdx = targetFileIdx;
      placeholderFt.featureIdx = featuresTargetRun.size();
      placeholderFt.intensity = 0.0;
      placeholderFt.rTime = *predRTimeIt;
      featuresTargetRun.push_back(placeholderFt);
      precursorLinks[queryIt->featureIdx] = 
          FeatureIdxMatch(placeholderFt.featureIdx, 1.0, placeholderFt.intensity);
      ++placeholdersAdded;
    }
  }
  
  if (Globals::VERB > 1) {
    std::cerr << "Added " << placeholdersAdded << " placeholders." << std::endl;
  }
} 

void FeatureAlignment::matchFeatures(
    const std::string& percolatorOutputFile,
    const DinosaurFeatureList& featuresQueryRun, 
    const DinosaurFeatureList& featuresTargetRun, 
    const std::vector<double>& predictedRTimesTargetRun, 
    float rTimeStdev) {  
  float rTimeTol = rTimeStdevTol_ * rTimeStdev;
  
  std::vector<std::string> percolatorArgs = percolatorArgs_;
  
  percolatorArgs.push_back("--results-psms");
  percolatorArgs.push_back(percolatorOutputFile);
  percolatorArgs.push_back("--decoy-results-psms");
  percolatorArgs.push_back(percolatorOutputFile + ".decoys");
  /*
  percolatorArgs.push_back("--tab-out");
  percolatorArgs.push_back(percolatorOutputFile + ".pin");
  */
  
  PercolatorAdapter percolatorAdapter;
  percolatorAdapter.parseOptions(percolatorArgs);
  percolatorAdapter.init(kLinkFeatureNames);
  
  DinosaurFeatureList::const_iterator queryIt, targetIt;
  std::vector<double>::const_iterator predRTimeIt;
  for (queryIt = featuresQueryRun.begin(), predRTimeIt = predictedRTimesTargetRun.begin(); 
       queryIt != featuresQueryRun.end(); ++queryIt, ++predRTimeIt) {
    std::vector<FeatureFeatureMatch> topTargets;
    
    DinosaurFeature queryFeature = *queryIt;
    queryFeature.rTime = *predRTimeIt;
    int label = 1;
    findMatches(label, rTimeTol, featuresTargetRun, queryFeature, topTargets);
    
    queryFeature.precMz += decoyOffset_;
    label = -1;
    findMatches(label, rTimeTol, featuresTargetRun, queryFeature, topTargets);
    
    std::sort(topTargets.begin(), topTargets.end());
    
    
    int numAdded = 0;
    std::vector<FeatureFeatureMatch>::const_iterator ffmIt;
    for (ffmIt = topTargets.begin(); 
          ffmIt != topTargets.end() && numAdded < maxFeatureCandidates_; 
          ++ffmIt, ++numAdded) {
      addLinkPsm(ffmIt->label, ffmIt->queryFeature, ffmIt->targetFeature, 
                 percolatorAdapter);
    }
  }
  
  percolatorAdapter.process();
}


void FeatureAlignment::addFeatureLinks(const std::string& percolatorOutputFile,
    DinosaurFeatureList& candidateFeaturesTargetRun,
    DinosaurFeatureList& featuresTargetRun,
    std::map<int, FeatureIdxMatch>& precursorLinks) {
  std::ifstream dataStream(percolatorOutputFile.c_str(), ios::in);
  
  std::cerr << "Links before "<< precursorLinks.size() << std::endl;
  
  int targetFileIdx = featuresTargetRun.begin()->fileIdx;
  std::string psmLine;
  getline(dataStream, psmLine); /* skip header */
  while (getline(dataStream, psmLine)) {
    TabReader reader(psmLine);
  
    std::string querySpecId = reader.readString();
    DinosaurFeature queryFt = parsePsmIdAsFeature(querySpecId);
    
    reader.readDouble(); /* score */
    double qval = reader.readDouble();
    double posteriorErrorProb = reader.readDouble();
    
    if (posteriorErrorProb > linkPEPThreshold_) break;
    
    if (posteriorErrorProb < precursorLinks[queryFt.featureIdx].posteriorErrorProb) {
      std::string peptide = reader.readString();
      std::string targetSpecId = peptide.substr(2, peptide.size() - 4);
      DinosaurFeature targetFt = parsePsmIdAsFeature(targetSpecId);
      if (targetFt.featureIdx < 0) {
        DinosaurFeature& candFt = candidateFeaturesTargetRun.at(
            getCandidatePos(targetFt.featureIdx));
        if (candFt.featureIdx < 0) {
          candFt.featureIdx = featuresTargetRun.size();
          candFt.fileIdx = targetFileIdx;
          featuresTargetRun.push_back(candFt);
        }
        targetFt = candFt;
      }
      precursorLinks[queryFt.featureIdx] = FeatureIdxMatch(targetFt.featureIdx, 
          posteriorErrorProb, targetFt.intensity);
    }
  }
  
  std::cerr << "Links after "<< precursorLinks.size() << std::endl;
}

void FeatureAlignment::findMatches(int label, double rTimeTol,
    const DinosaurFeatureList& featuresTargetRun, 
    const DinosaurFeature& queryFeature, 
    std::vector<FeatureFeatureMatch>& topTargets) {
  DinosaurFeatureList::const_iterator lowerBound = 
      featuresTargetRun.getPrecMzIterator(queryFeature.precMz*(1 - ppmTol_*1e-6));
  DinosaurFeatureList::const_iterator upperBound = 
      featuresTargetRun.getPrecMzIterator(queryFeature.precMz*(1 + ppmTol_*1e-6));
  DinosaurFeatureList::const_iterator targetIt;
  for (targetIt = lowerBound; targetIt != upperBound; ++targetIt) {
    float rTimeDiff = std::abs(targetIt->rTime - queryFeature.rTime);
    if (targetIt->charge == queryFeature.charge && rTimeDiff <= rTimeTol) {
      DinosaurFeature targetFeature = *targetIt;
      float ppmDiff = getPpmDiff(targetIt->precMz, queryFeature.precMz);
      float score = rTimeDiff / rTimeTol + ppmDiff / ppmTol_; /* score <= 1.0 */
      /* put placeholder features in the back of the queue */
      if (targetIt->intensity == 0.0) score += 1.0; 
      
      topTargets.push_back(FeatureFeatureMatch(score, queryFeature, targetFeature, label));
    }
  }
}
    
void FeatureAlignment::addLinkPsm(int label,
    const DinosaurFeature& queryFeature, 
    const DinosaurFeature& targetFeature, 
    PercolatorAdapter& percolatorAdapter) {
  int scannr = queryFeature.featureIdx;
  std::string querySpecId = convertFeatureToPsmId(queryFeature);
  std::string targetSpecId = convertFeatureToPsmId(targetFeature);
  std::string peptide = "A." + targetSpecId + ".A";
  
  double ppmDiff = getPpmDiff(targetFeature.precMz, queryFeature.precMz);
  double rTimeDiff = std::abs(queryFeature.rTime - targetFeature.rTime);
  /* double xicRateDiff = getXicRateDiff(targetFeature, queryFeature); */
  int queryIsPlaceHolder = (queryFeature.intensity == 0.0) ? 1 : 0;
  int targetIsPlaceHolder = (targetFeature.intensity == 0.0) ? 1 : 0;
  
  std::vector<double> features;
  features.push_back(ppmDiff);
  features.push_back(rTimeDiff);
  features.push_back(targetFeature.precMz);
  features.push_back(targetFeature.rTime);
  /* features.push_back(xicRateDiff); */
  features.push_back(queryIsPlaceHolder);
  features.push_back(targetIsPlaceHolder);
  for (int charge = 1; charge <= 4; ++charge) {
    int isCharge = (queryFeature.charge == charge || 
                    (queryFeature.charge > 4 && charge == 4)) ? 1 : 0;
    features.push_back(isCharge);
  }
  
  percolatorAdapter.addPsm(querySpecId, label, scannr, peptide, features);
}

double FeatureAlignment::getXicRateDiff(const DinosaurFeature& queryFt, 
    const DinosaurFeature& targetFt) {
  double queryIntensity = (queryFt.intensity == 0.0) ? targetFt.intensity : queryFt.intensity;
  double targetIntensity = (targetFt.intensity == 0.0) ? queryFt.intensity : targetFt.intensity;
  if (queryIntensity * targetIntensity == 0.0) {
    queryIntensity = targetIntensity = 1e6;
  }
  double targetXicRate = getXicRate(targetIntensity, targetFt.rtStart, targetFt.rtEnd);
  double queryXicRate = getXicRate(queryIntensity, queryFt.rtStart, queryFt.rtEnd);
  return std::abs(queryXicRate - targetXicRate);
}

double FeatureAlignment::getXicRate(const double intensity, 
    const double rtStart, const double rtEnd) {
  return log(intensity) / ((rtEnd - rtStart)*60.0);
}

void FeatureAlignment::mbrMatchFeatures(
    const std::map<int, FeatureIdxMatch>& precursorLinks,
    const std::string& percolatorOutputFile,
    const std::string& targetMzMLFile,
    const DinosaurFeatureList& featuresQueryRun, 
    const DinosaurFeatureList& featuresTargetRun, 
    const std::vector<double>& predictedRTimesTargetRun, 
    float rTimeStdev,
    DinosaurFeatureList& candidateFeaturesTargetRun) {
  float rTimeTol = rTimeStdevTol_ * rTimeStdev;
  
  /* TODO: put these in a separate folder */
  std::string targetFile = percolatorOutputFile + ".dinosaur_targets.tsv";
  if (Globals::fileIsEmpty(targetFile)) {
    std::ofstream targetFileStream(targetFile.c_str(), ios::out);
    targetFileStream << "mz\tcharge\tmzDiff\trtStart\trtEnd\tminApexInt\tid" << std::endl;
    
    DinosaurFeatureList::const_iterator queryIt;
    std::vector<double>::const_iterator predRTimeIt;
    int targetFileIdx = featuresTargetRun.begin()->fileIdx;
    size_t placeholdersAdded = 0u;
    for (queryIt = featuresQueryRun.begin(), predRTimeIt = predictedRTimesTargetRun.begin(); 
         queryIt != featuresQueryRun.end(); ++queryIt, ++predRTimeIt) {
      if (precursorLinks.find(queryIt->featureIdx) == precursorLinks.end() || 
          precursorLinks.find(queryIt->featureIdx)->second.posteriorErrorProb >= linkPEPMbrSearchThreshold_) {
        std::ostringstream commonStream;
        commonStream << queryIt->charge << "\t" <<
                        queryIt->precMz * ppmTol_ * 1e-6 << "\t" <<
                        *predRTimeIt - rTimeTol << "\t" <<
                        *predRTimeIt + rTimeTol << "\t" <<
                        10000 << "\t" <<
                        queryIt->featureIdx << std::endl;
        
        targetFileStream << queryIt->precMz << "\t" << commonStream.str() << 
            queryIt->precMz + decoyOffset_ << "\t" << commonStream.str();
      }
    }
    targetFileStream.close();
  }
  
  boost::filesystem::path targetMzMLFilePath(targetMzMLFile);
  std::string outputDir = percolatorOutputFile + "_dinosaur"; /* TODO: put these in a separate folder */
  std::string dinosaurOutputFile = outputDir + "/" + 
      targetMzMLFilePath.stem().string() + ".targets.csv";
  if (Globals::fileIsEmpty(dinosaurOutputFile)) {
    boost::filesystem::path dinosaurOutputDir(outputDir);
    boost::system::error_code returnedError;
    boost::filesystem::create_directories(dinosaurOutputDir, returnedError);
    if (!boost::filesystem::exists(dinosaurOutputDir)) {
      std::ostringstream oss;
      oss << "Error: could not create output directory at " << 
             dinosaurOutputDir.string() << std::endl;
      throw MyException(oss.str());
    }
    
    int rc = DinosaurIO::runDinosaurTargeted(outputDir, targetMzMLFile, targetFile);
    if (rc != EXIT_SUCCESS) {
      std::ostringstream oss;
      oss << "Dinosaur failed with exit code " << rc << 
             ". Terminating.." << std::endl;
      throw MyException(oss.str());
    }
  }
  
  processDinosaurTargets(precursorLinks, dinosaurOutputFile, 
      percolatorOutputFile, featuresQueryRun, featuresTargetRun, 
      predictedRTimesTargetRun, candidateFeaturesTargetRun);
}
    
void FeatureAlignment::processDinosaurTargets(
    const std::map<int, FeatureIdxMatch>& precursorLinks,
    const std::string& dinosaurTargetOutputFile, 
    const std::string& percolatorOutputFile, 
    const DinosaurFeatureList& featuresQueryRun,
    const DinosaurFeatureList& featuresTargetRun,
    const std::vector<double>& predictedRTimesTargetRun,
    DinosaurFeatureList& candidateFeaturesTargetRun) {
  std::vector<std::string> percolatorArgs = percolatorArgs_;
  
  PercolatorAdapter percolatorAdapter;
  bool runPercolator = false;
  if (Globals::fileIsEmpty(percolatorOutputFile)) {
    runPercolator = true;
    percolatorArgs.push_back("--results-psms");
    percolatorArgs.push_back(percolatorOutputFile);
    percolatorArgs.push_back("--decoy-results-psms");
    percolatorArgs.push_back(percolatorOutputFile + ".decoys");
    /*
    percolatorArgs.push_back("--tab-out");
    percolatorArgs.push_back(percolatorOutputFile + ".pin");
    */
    
    percolatorAdapter.parseOptions(percolatorArgs);  
    percolatorAdapter.init(kLinkFeatureNames);
  }
  
  std::ifstream dinosaurOutputStream(dinosaurTargetOutputFile.c_str(), ios::in);
  
  std::string targetLine;
  getline(dinosaurOutputStream, targetLine); /* skip header */
  
  DinosaurFeatureList::const_iterator queryIt;
  std::vector<double>::const_iterator predRTimeIt;
  for (queryIt = featuresQueryRun.begin(), predRTimeIt = predictedRTimesTargetRun.begin(); 
       queryIt != featuresQueryRun.end(); ++queryIt, ++predRTimeIt) {
    if (precursorLinks.find(queryIt->featureIdx) == precursorLinks.end() || 
        precursorLinks.find(queryIt->featureIdx)->second.posteriorErrorProb >= linkPEPMbrSearchThreshold_) {
      getline(dinosaurOutputStream, targetLine);
      DinosaurFeature targetFeature = DinosaurIO::parseDinosaurFeatureRow(targetLine);
      targetFeature.featureIdx = featuresTargetRun.getFeatureIdx(targetFeature);
      targetFeature.fileIdx = -1;
      if (targetFeature.intensity > 0.0 && targetFeature.featureIdx == -1) {
        targetFeature.featureIdx = candidateFeaturesTargetRun.getFeatureIdx(targetFeature);
        if (targetFeature.featureIdx == -1) {
          targetFeature.featureIdx = getCandidateIdx(candidateFeaturesTargetRun.size());
          candidateFeaturesTargetRun.push_back(targetFeature);
        }
      }
      
      getline(dinosaurOutputStream, targetLine);
      if (runPercolator) {
        DinosaurFeature decoyFeature = DinosaurIO::parseDinosaurFeatureRow(targetLine);
        decoyFeature.featureIdx = -1;
        decoyFeature.fileIdx = -1;
        
        DinosaurFeature queryFeature = *queryIt;
        queryFeature.rTime = *predRTimeIt;
        if (targetFeature.intensity > 0.0)
          addLinkPsm(1, queryFeature, targetFeature, percolatorAdapter);
        
        queryFeature.precMz += decoyOffset_;
        if (decoyFeature.intensity > 0.0)
          addLinkPsm(-1, queryFeature, decoyFeature, percolatorAdapter);
      }
    }
  }
  
  if (runPercolator) {    
    percolatorAdapter.process();
  }
}

std::string FeatureAlignment::convertFeatureToPsmId(const DinosaurFeature& feature) {
  return boost::lexical_cast<std::string>(feature.featureIdx) + "_" + 
         boost::lexical_cast<std::string>(feature.precMz) + "_" + 
         boost::lexical_cast<std::string>(feature.rTime) + "_" + 
         boost::lexical_cast<std::string>(feature.intensity) + "_" + 
         boost::lexical_cast<std::string>(feature.charge);
}

DinosaurFeature FeatureAlignment::parsePsmIdAsFeature(const std::string& specId) {
  DinosaurFeature ft;
  
  size_t prev = 0u;
  size_t next = specId.find_first_of("_", prev);
  ft.featureIdx = boost::lexical_cast<int>(specId.substr(prev, next - prev));
  prev = next+1;
  
  next = specId.find_first_of("_", prev);
  ft.precMz = boost::lexical_cast<float>(specId.substr(prev, next - prev));
  prev = next+1;
  
  next = specId.find_first_of("_", prev);
  ft.rTime = boost::lexical_cast<float>(specId.substr(prev, next - prev));
  prev = next+1;
  
  next = specId.find_first_of("_", prev);
  ft.intensity = boost::lexical_cast<float>(specId.substr(prev, next - prev));
  prev = next+1;
  
  next = specId.find_first_of("_", prev);
  ft.charge = boost::lexical_cast<int>(specId.substr(prev, next - prev));
  
  return ft;
}

} /* namespace quandenser */
