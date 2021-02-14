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

#include "AlignRetention.h"

namespace quandenser {

bool operator==(const FilePair& l, const FilePair& r) {
  return l.fileIdx1 == r.fileIdx1 && l.fileIdx2 == r.fileIdx2;
}

bool operator<(const FilePair& l, const FilePair& r) {
  return (l.fileIdx1 < r.fileIdx1 || (l.fileIdx1 == r.fileIdx1 && l.fileIdx2 < r.fileIdx2));
}

bool operator==(const RTimePair& l, const RTimePair& r) {
  return l.rTime1 == r.rTime1 && l.rTime2 == r.rTime2;
}

bool operator<(const RTimePair& l, const RTimePair& r) {
  return (l.rTime1 < r.rTime1 || (l.rTime1 == r.rTime1 && l.rTime2 < r.rTime2));
}

void AlignRetention::getAlignModelsTree() {
  getAlignModels();

  int maxFileIdx = 0;
  RTimePairs::iterator filePairIt;
  std::vector<std::pair<float, FilePair> > sortedWeights;
  for (filePairIt = rTimePairs_.begin(); filePairIt != rTimePairs_.end(); ++filePairIt) {
    if (rmses_.find(filePairIt->first) == rmses_.end()) continue;

    FilePair revFilePair = filePairIt->first.getRevFilePair();
    double rmse1 = rmses_[filePairIt->first];
    double rmse2 = rmses_[revFilePair];

    maxFileIdx = std::max(maxFileIdx, revFilePair.fileIdx1);
    maxFileIdx = std::max(maxFileIdx, revFilePair.fileIdx2);

    float rmseComb = std::sqrt(rmse1*rmse1 + rmse2*rmse2);

    if (!isinf(rmseComb) && !isnan(rmseComb) && rmse1 > 0 && rmse2 > 0) {
      sortedWeights.push_back(std::make_pair(rmseComb, filePairIt->first));
    }
  }
  int numFiles = maxFileIdx + 1;

  std::sort(sortedWeights.begin(), sortedWeights.end());

  for (int i = 0; i < numFiles; ++i) {
    fileGraphNodes_.push_back(FileGraphNode(i));
  }

  std::set<int> addedNodes;
  std::set<FilePair> addedLinks;
  if (sortedWeights.size() > 0) {
    addedNodes.insert(sortedWeights.front().second.fileIdx1);
  }
  while (addedNodes.size() < numFiles && sortedWeights.size() > 0) {
    bool inserted = false;
    for (std::vector<std::pair<float, FilePair> >::const_iterator it = sortedWeights.begin();
          it != sortedWeights.end(); ++it) {
      int fileIdx1 = it->second.fileIdx1;
      int fileIdx2 = it->second.fileIdx2;
      /* exactly 1 of the fileIdxs has to be in the addedNodes set */
      if ((addedNodes.find(fileIdx1) != addedNodes.end()) != (addedNodes.find(fileIdx2) != addedNodes.end())) {
        inserted = true;
        addedNodes.insert(fileIdx1);
        addedNodes.insert(fileIdx2);
        fileGraphNodes_[fileIdx1].addNeighbor(fileIdx2);
        fileGraphNodes_[fileIdx2].addNeighbor(fileIdx1);

        addedLinks.insert(it->second);

        if (Globals::VERB > 2) {
          std::cerr << "Inserting link " << fileIdx1 << " to "
                    << fileIdx2 << " with rmse " << it->first << std::endl;
        }

        break;
      }
    }
    if (!inserted) break;
  }

  if (addedNodes.size() != numFiles) {
    std::ostringstream oss;
    oss << "ERROR: Could not create minimum spanning tree for alignments." << std::endl
        << "  Not enough overlap of MS2 clusters for the following runs:";
    for (int i = 0; i < numFiles; ++i) {
      if (addedNodes.find(i) == addedNodes.end()) {
        oss << " " << i;
      }
    }
    oss << std::endl << "  Terminating.." << std::endl;
    throw MyException(oss.str());
  }

  std::set<FilePair>::iterator addedLinkIt;
  for (addedLinkIt = addedLinks.begin(); addedLinkIt != addedLinks.end(); ++addedLinkIt) {
    FilePair filePair = *addedLinkIt;
    std::vector<RTimePair> rTimePairs = rTimePairs_[filePair];
    
    std::vector<double> medianRTimesRun1, medianRTimesRun2;
    getKnots(rTimePairs, medianRTimesRun1, medianRTimesRun2);
    
    alignments_[filePair].setData(medianRTimesRun1, medianRTimesRun2);
    alignments_[filePair].roughnessPenaltyIRLS();
    bool reversedPair = false;
    float rmse1 = getRMSE(alignments_[filePair], rTimePairs, reversedPair);
    alignments_[filePair].setRmse(rmse1);

    FilePair revFilePair = filePair.getRevFilePair();
    alignments_[revFilePair].setData(medianRTimesRun2, medianRTimesRun1);
    alignments_[revFilePair].roughnessPenaltyIRLS();
    reversedPair = true;
    float rmse2 = getRMSE(alignments_[revFilePair], rTimePairs, reversedPair);
    alignments_[revFilePair].setRmse(rmse2);
  }
  
  rTimePairs_.clear();
  rmses_.clear();
}

void AlignRetention::getAlignModels() {  
  /* "touch" all keys to allow OMP to concurrently fill it later */
  RTimePairs::iterator filePairItReserve;
  for (filePairItReserve = rTimePairs_.begin(); filePairItReserve != rTimePairs_.end(); ++filePairItReserve) {
    rmses_[filePairItReserve->first];
    FilePair revFilePair = filePairItReserve->first.getRevFilePair();
    rmses_[revFilePair];
  }

#pragma omp parallel for schedule(dynamic, 100)
  for (int i = 0; i < rTimePairs_.size(); ++i) {
    RTimePairs::iterator filePairIt = rTimePairs_.begin();
    std::advance(filePairIt, i);
    
    std::vector<double> medianRTimesRun1, medianRTimesRun2;

    if (filePairIt->second.size() < 250) {
      if (Globals::VERB > 2) {
        std::cerr << "Skipping aligning runs: " << filePairIt->first.fileIdx1 << " " << filePairIt->first.fileIdx2
            << ". Fewer than 250 retention time pairs: " << filePairIt->second.size() << std::endl;
      }
      continue;
    }

    getKnots(filePairIt->second, medianRTimesRun1, medianRTimesRun2);
    
    if (medianRTimesRun1.size() < 50) {
      if (Globals::VERB > 2) {
        std::cerr << "Skipping aligning runs: " << filePairIt->first.fileIdx1 << " " << filePairIt->first.fileIdx2
            << ". Fewer than 50 retention time knots: " << medianRTimesRun1.size() <<  std::endl;
      }
      continue;
    }
    
    SplineRegression splineRegression;
    splineRegression.setData(medianRTimesRun1, medianRTimesRun2);
    splineRegression.roughnessPenaltyIRLS();
    bool reversedPair = false;
    float rmse1 = getRMSE(splineRegression, filePairIt->second, reversedPair);
    rmses_[filePairIt->first] = rmse1;
    
    SplineRegression splineRegressionReverse;
    FilePair revFilePair = filePairIt->first.getRevFilePair();
    splineRegressionReverse.setData(medianRTimesRun2, medianRTimesRun1);
    splineRegressionReverse.roughnessPenaltyIRLS();
    reversedPair = true;
    float rmse2 = getRMSE(splineRegressionReverse, filePairIt->second, reversedPair);
    rmses_[revFilePair] = rmse2;

    if (Globals::VERB > 2) {
      float rmseComb = std::sqrt(rmse1*rmse1 + rmse2*rmse2);
      std::cerr << "Aligned runs: " << filePairIt->first.fileIdx1 << " " << filePairIt->first.fileIdx2
          << ": rmseComb = " << rmseComb <<  " rmse1 = " << rmse1
          << " rmse2 = " << rmse2 << " numPairs = " << filePairIt->second.size() << std::endl;
    }
  }
}

/* TODO: use cross validation to select an appropriate number of bins */
/* TODO: use pair with median error in a certain range around the medianIdx. This will be an actual datapoint, but somewhat protected against outliers */
void AlignRetention::getKnots(std::vector<RTimePair>& rTimePairsSingleFilePair, std::vector<double>& medianRTimesRun1, std::vector<double>& medianRTimesRun2) {
  std::sort(rTimePairsSingleFilePair.begin(), rTimePairsSingleFilePair.end());
  rTimePairsSingleFilePair.erase(std::unique(rTimePairsSingleFilePair.begin(), rTimePairsSingleFilePair.end() ), rTimePairsSingleFilePair.end() );

  std::vector<double> rTimesRun1, rTimesRun2;
  std::vector<RTimePair>::const_iterator rtPair;
  for (rtPair = rTimePairsSingleFilePair.begin(); rtPair != rTimePairsSingleFilePair.end(); ++rtPair) {
    rTimesRun1.push_back(rtPair->rTime1);
    rTimesRun2.push_back(rtPair->rTime2);
  }
  std::sort(rTimesRun1.begin(), rTimesRun1.end());
  std::sort(rTimesRun2.begin(), rTimesRun2.end());
  
  int numBins = 100;
  float binSize = static_cast<float>(rTimePairsSingleFilePair.size()) / numBins;
  for (int bin = 0; bin < numBins; ++bin) {
    size_t medianIdx = static_cast<int>(std::floor(binSize * (bin + 0.5f)));
    //RTimePair rtPair = rTimePairsSingleFilePair.at(medianIdx);
    if (medianIdx < rTimesRun1.size()
        && (medianRTimesRun1.size() == 0
          || (medianRTimesRun1.back() != rTimesRun1.at(medianIdx) && medianRTimesRun2.back() != rTimesRun2.at(medianIdx) ) ) ) {
      medianRTimesRun1.push_back(rTimesRun1.at(medianIdx));
      medianRTimesRun2.push_back(rTimesRun2.at(medianIdx));
    }
  }
}

/*
  Robust regression: MAD*1.4826 (https://en.wikipedia.org/wiki/Robust_measures_of_scale)
  Alternatively, we could try this: Huber (http://www.statsmodels.org/dev/generated/statsmodels.robust.scale.Huber.html)
 */
float AlignRetention::getRMSE(SplineRegression& alignment, std::vector<RTimePair>& rTimePairsSingleFilePair, bool reversedPair) {
  std::vector<double> rTimesRun1, rTimesRun2, predictedRTimesRun2;
  std::vector<RTimePair>::const_iterator rtPair;
  for (rtPair = rTimePairsSingleFilePair.begin(); rtPair != rTimePairsSingleFilePair.end(); ++rtPair) {
    rTimesRun1.push_back(rtPair->rTime1);
    rTimesRun2.push_back(rtPair->rTime2);
  }
  if (reversedPair) rTimesRun1.swap(rTimesRun2);
  alignment.predict(rTimesRun1, predictedRTimesRun2);

  std::vector<double> absErrors;
  std::vector<double>::const_iterator run2RtIt, predRtIt;
  for (run2RtIt = rTimesRun2.begin(), predRtIt = predictedRTimesRun2.begin(); run2RtIt != rTimesRun2.end(); ++run2RtIt, ++predRtIt) {
    absErrors.push_back(std::abs(*run2RtIt - *predRtIt));
    //std::cerr << std::abs(*run2RtIt - *predRtIt) << " " << *run2RtIt << " " << *predRtIt << std::endl;
  }
  /*
  size_t numBins = 100;
  float binSize = float(absErrors.size()) /  numBins;
  for (int i = 0; i < numBins; ++i) {
    std::sort(absErrors.begin() + size_t(i * binSize), absErrors.begin() + std::min(absErrors.size(), size_t((i+1) * binSize)));
    std::cerr << i << " " << absErrors.at(size_t(i * binSize) + size_t(binSize) / 2) * 1.4826 << std::endl;
  }
  */
  std::sort(absErrors.begin(), absErrors.end());
  float rmse = absErrors.at(absErrors.size() / 2) * 1.4826;
  return rmse;
}

void AlignRetention::createMinDepthTree(std::vector<std::pair<int, FilePair> >& featureAlignmentQueue) {
  int furthestNode = breadthFirstSearch(0);
  furthestNode = breadthFirstSearch(furthestNode);

  int rootNode = getMinimumDepthRoot(furthestNode);
  breadthFirstSearch(rootNode);

  /* we can process all feature alignments of the same round in parallel, before moving onto the next round */
  for (std::vector<FileGraphNode>::const_iterator nodeIt = fileGraphNodes_.begin(); nodeIt != fileGraphNodes_.end(); ++nodeIt) {
    int nodeDepth = nodeIt->getDepth();
    int round = maxDepth_ - nodeDepth;

    if (nodeDepth > 0) {
      FilePair filePair(nodeIt->getFileIdx(), nodeIt->getParent());
      if (Globals::VERB > 2) {
        std::cerr << "Inserting feature alignment " << filePair.fileIdx1 
          << "->" << filePair.fileIdx2 << " to the queue for round " 
          << round << "." << std::endl;
      }
      featureAlignmentQueue.push_back(std::make_pair(round, filePair));
    }

    int childRound = maxDepth_ + nodeDepth;
    for (std::set<int>::const_iterator childNodeIt = nodeIt->getChildrenItBegin(); childNodeIt != nodeIt->getChildrenItEnd(); ++childNodeIt) {
      FilePair childFilePair(nodeIt->getFileIdx(), *childNodeIt);
      if (Globals::VERB > 2) {
        std::cerr << "Inserting feature alignment " << childFilePair.fileIdx1 
          << "->" << childFilePair.fileIdx2 << " to the queue for round " 
          << childRound << "." << std::endl;
      }
      featureAlignmentQueue.push_back(std::make_pair(childRound, childFilePair));
    }
  }
  std::sort(featureAlignmentQueue.begin(), featureAlignmentQueue.end());
}

int AlignRetention::breadthFirstSearch(int startNode) {
  /* contains pairs (-1*insertion_depth, node_to_investigate). The -1 is
     necessary because the default priority queue takes high values first */
  std::priority_queue<std::pair<int, int> > searchQueue;

  int depth = 0;
  int parent = -1;
  std::set<int> newChildren = fileGraphNodes_[startNode].getNeighbors();
  fileGraphNodes_[startNode].insertInTree(parent, newChildren, depth);
  searchQueue.push(std::make_pair(-1*(depth + 1), startNode));

  int lastNode = -1;
  while (searchQueue.size() > 0) {
    std::pair<int, int> depthNodePair = searchQueue.top();
    searchQueue.pop();
    depth = -1*depthNodePair.first;
    startNode = depthNodePair.second;

    std::set<int>::const_iterator childIt;
    for (childIt = fileGraphNodes_[startNode].getChildrenItBegin(); childIt != fileGraphNodes_[startNode].getChildrenItEnd(); ++childIt) {
      newChildren = fileGraphNodes_[*childIt].getNeighbors();
      newChildren.erase(startNode);

      fileGraphNodes_[*childIt].insertInTree(startNode, newChildren, depth);
      searchQueue.push(std::make_pair(-1*(depth + 1), *childIt));
    }
    lastNode = startNode;
  }
  return lastNode;
}

int AlignRetention::getMinimumDepthRoot(int lastNode) {
  maxDepth_ = static_cast<int>(std::ceil((fileGraphNodes_[lastNode].getDepth()+1) / 2.0f));
  while (fileGraphNodes_[lastNode].getDepth() > maxDepth_) {
    lastNode = fileGraphNodes_[lastNode].getParent();
  }
  return lastNode;
}

int AlignRetention::saveState(const std::string& alignFilePath){
  std::ofstream outfile(alignFilePath.c_str());
  /* To save:
  RTimePairs rTimePairs_;
  std::map<FilePair, SplineRegression> alignments_;
  std::vector<FileGraphNode> fileGraphNodes_;
  int maxDepth_;
  */

  // alignment_
  if (outfile.is_open()) {
    std::map<FilePair, SplineRegression>::iterator alignmentIt;
    for (alignmentIt = alignments_.begin(); alignmentIt != alignments_.end(); alignmentIt++) {

      // Save filepair + rmse
      FilePair filePair = alignmentIt->first;
      SplineRegression spline = alignmentIt->second;
      outfile << filePair.fileIdx1 << '\t' << filePair.fileIdx2 << '\t' << spline.getRmse() << '\t';

      // Save data
      std::pair<std::vector<double>, std::vector<double> > vector_pair = spline.getData();
      std::vector<double> x = vector_pair.first;
      std::vector<double> y = vector_pair.second;
      
      outfile << x.size() << '\t';
      for (int i = 0; i < x.size(); i++) {
        outfile << x[i] << '\t';
      }
      outfile << y.size() << '\t';
      for (int i = 0; i < y.size(); i++) {
        outfile << y[i] << '\t';
      }

      outfile << std::endl;
    }
    outfile.close();
    return EXIT_SUCCESS;
  } else {
    std::cerr << "ERROR: Could not open " << alignFilePath
        << " for writing the retention time alignments." << std::endl;
    return EXIT_FAILURE;
  }
}

void AlignRetention::loadState(const std::string& alignFilePath, const FilePair& alignedFilePair) {
  // Load state of featureAlignmentQueue
  std::cout << "Loading state of alignment" << std::endl;
  std::ifstream infile(alignFilePath.c_str());
  SplineRegression spline;
  int fileidx1;
  int fileidx2;
  FilePair dummyFilePair(-1, -1);
  while (infile >> fileidx1 >> fileidx2) {
    FilePair filePair(fileidx1, fileidx2);
    
    // Read rmse
    float rmse;
    infile >> rmse;

    // Read x vector
    unsigned long int x_size;
    infile >> x_size;
    std::vector<double> x;
    for (unsigned long int i = 0; i < x_size; ++i) {
      std::string d_string;  // To string, because nan could be in file
      infile >> d_string;
      double d = atof(d_string.c_str());
      x.push_back(d);
    }

    // Read y vector
    unsigned long int y_size;
    infile >> y_size;
    std::vector<double> y;
    for (unsigned long int i = 0; i < y_size; ++i) {
      std::string d_string;  // To string, because nan could be in file
      infile >> d_string;
      double d = atof(d_string.c_str());
      y.push_back(d);
    }
    
    if (!(alignedFilePair == dummyFilePair || filePair == alignedFilePair)) continue;
    
    spline.setData(x,y);
    spline.roughnessPenaltyIRLS();
    spline.predict(x,y);
    spline.setRmse(rmse);
    alignments_[filePair] = spline;
  }

  /*
    // rTimePairs_
    std::vector<double> medianRTimesRun1, medianRTimesRun2;
    getKnots(filePairIt->second, medianRTimesRun1, medianRTimesRun2);
    alignments_[filePairIt->first].setData(medianRTimesRun1, medianRTimesRun2);
    alignments_[filePairIt->first].roughnessPenaltyIRLS();
    bool reversedPair = false;
    float rmse1 = getRMSE(alignments_[filePairIt->first], filePairIt->second, reversedPair);
    alignments_[filePairIt->first].setRmse(rmse1);
    FilePair revFilePair = filePairIt->first.getRevFilePair();
    alignments_[revFilePair].setData(medianRTimesRun2, medianRTimesRun1);
    alignments_[revFilePair].roughnessPenaltyIRLS();
    reversedPair = true;
    float rmse2 = getRMSE(alignments_[revFilePair], filePairIt->second, reversedPair);
    alignments_[revFilePair].setRmse(rmse2);
  */

  std::cout << "State loaded" << std::endl;
  infile.close();
}

} /* namespace quandenser */
