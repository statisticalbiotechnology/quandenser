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

#ifndef QUANDENSER_ALIGNRETENTION_H_
#define QUANDENSER_ALIGNRETENTION_H_

#include <vector>
#include <map>
#include <set>
#include <queue>

#include "SplineRegression.h"
#include "MyException.h"

namespace quandenser {

struct FilePair {
  FilePair(int f1, int f2) : fileIdx1(f1), fileIdx2(f2) {}
  int fileIdx1, fileIdx2;

  FilePair getRevFilePair() const { return FilePair(fileIdx2, fileIdx1); }
};

bool operator==(const FilePair& l, const FilePair& r);
bool operator<(const FilePair& l, const FilePair& r);

struct RTimePair {
  RTimePair(float rt1, float rt2) : rTime1(rt1), rTime2(rt2) {}
  float rTime1, rTime2;
};

bool operator==(const RTimePair& l, const RTimePair& r);
bool operator<(const RTimePair& l, const RTimePair& r);

typedef std::map<FilePair, std::vector<RTimePair> > RTimePairs;

class FileGraphNode {
 public:
  FileGraphNode(int fileIdx) : fileIdx_(fileIdx), depth_(-1), parentFileIdx_(-1),
    children_(), neighbors_() {}

  inline void addNeighbor(int idx) { neighbors_.insert(idx); }
  inline void insertInTree(int parentFileIdx, std::set<int>& children, int depth) {
    parentFileIdx_ = parentFileIdx;
    children_ = children;
    depth_ = depth;
  }

  inline int getDepth() const { return depth_; }
  inline int getParent() const { return parentFileIdx_; }
  inline int getFileIdx() const { return fileIdx_; }
  inline std::set<int> getNeighbors() { return neighbors_; }

  inline std::set<int>::const_iterator getChildrenItBegin() const { return children_.begin(); }
  inline std::set<int>::const_iterator getChildrenItEnd() const { return children_.end(); }
 protected:
  int fileIdx_;
  int depth_;
  int parentFileIdx_;
  std::set<int> children_;
  std::set<int> neighbors_;
};

class AlignRetention {
 public:
  AlignRetention() : rTimePairs_() {}

  inline RTimePairs& getRTimePairsRef() { return rTimePairs_; }
  inline int getMaxDepth() const { return maxDepth_; }

  void getAlignModelsTree();
  void createMinDepthTree(std::vector<std::pair<int, FilePair> >& featureAlignmentQueue);
  int saveState(const std::string& alignFilePath);
  void loadState(const std::string& alignFilePath);

  SplineRegression& getAlignment(FilePair& filePair) { return alignments_[filePair]; }

 protected:
  RTimePairs rTimePairs_;
  std::map<FilePair, SplineRegression> alignments_;
  std::vector<FileGraphNode> fileGraphNodes_;
  int maxDepth_;

  void getAlignModels();
  void getKnots(std::vector<RTimePair>& rTimePairsSingleFilePair,
      std::vector<double>& medianRTimesRun1,
      std::vector<double>& medianRTimesRun2);
  float getRMSE(SplineRegression& alignment,
      std::vector<RTimePair>& rTimePairsSingleFilePair, bool reversedPair);

  int breadthFirstSearch(int startNode);
  int getMinimumDepthRoot(int lastNode);

};

} /* namespace quandenser */

#endif /* QUANDENSER_ALIGNRETENTION_H_ */
