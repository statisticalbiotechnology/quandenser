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

namespace quandenser {

/* Dijkstra's algorithm */
template <class Node>
void SimilarityMatrix<Node>::computeShortestPathsFromSource(Node key1, 
    std::map<Node, float>& similarities) {
  
  std::priority_queue<std::pair<float, Node> > searchQueue;
  searchQueue.push(std::make_pair(1.0, key1));
  
  while (searchQueue.size() > 0) {
    std::pair<float, Node> simNodePair = searchQueue.top();
    float sim = simNodePair.first;
    Node node1 = simNodePair.second;
    searchQueue.pop();
    if (similarities.find(node1) == similarities.end()) {
      similarities[node1] = sim;
      if (similarities.size() == matrix_.size()) break;
      
      typename std::map<Node, float>::const_iterator node2It;
      std::map<Node, float>& row = matrix_[node1];
      for (node2It = row.begin(); node2It != row.end(); ++node2It) {
        if (similarities.find(node2It->first) == similarities.end()) {
          searchQueue.push(std::make_pair(sim * node2It->second, node2It->first));
        }
      }
    }
  }
}

template <class Node>
void SimilarityMatrix<Node>::updateLink(Node key1, Node key2, float sim) {
  if (sim > getSim(key1, key2)) {
    matrix_[key1][key2] = sim;
    matrix_[key2][key1] = sim;
  }
}

} /* namespace quandenser */
