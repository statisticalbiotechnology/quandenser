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

#ifndef QUANDENSER_SIMILARITYMATRIX_H_
#define QUANDENSER_SIMILARITYMATRIX_H_

#include <vector>
#include <map>
#include <utility>
#include <queue>

namespace quandenser {

template <class Node> class SimilarityMatrix {
 protected:
  typedef typename std::map<Node, std::map<Node, float> >::iterator iterator;
  typedef typename std::map<Node, std::map<Node, float> >::const_iterator const_iterator;
  
 public:
  SimilarityMatrix(){};
  ~SimilarityMatrix(){};
  
  void computeShortestPathsFromSource(
    std::map<Node, float>& similarities);
  void updateLink(Node key1, Node key2, float sim);
  inline float getSim(Node key1, Node key2) {
    return (key1 == key2) ? 1.0 : matrix_[key1][key2];
  }
  
  /* delegation functions */
  inline iterator begin() { return matrix_.begin(); }
  inline iterator end() { return matrix_.end(); }
  inline const_iterator begin() const { return matrix_.begin(); }
  inline const_iterator end() const { return matrix_.end(); }
  
 protected:
  std::map<Node, std::map<Node, float> > matrix_;
  
};

} /* namespace quandenser */

#include "SimilarityMatrix.tpp"

#endif /* QUANDENSER_SIMILARITYMATRIX_H_ */
