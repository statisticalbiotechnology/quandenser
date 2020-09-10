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

#ifndef QUANDENSER_DINOSAURFEATURELIST_H_
#define QUANDENSER_DINOSAURFEATURELIST_H_

#include <vector>
#include <sstream>

#include <boost/unordered/unordered_map.hpp>

#include "maracluster/src/BinaryInterface.h"

#include "Globals.h"
#include "DinosaurFeature.h"

namespace quandenser {

class DinosaurFeatureList {
 protected:
  typedef std::vector<DinosaurFeature> FeatureList;

  struct DinosaurFeatureLessPrecMzLower {
    bool operator ()(const DinosaurFeature& ms, const float f) const {
      return ms.precMz < f;
    }
  };

  struct DinosaurFeatureEqual : std::binary_function<DinosaurFeature, DinosaurFeature, bool> {
    bool operator()(const DinosaurFeature& x, const DinosaurFeature& y) const {
      return x.precMz == y.precMz && x.rTime == y.rTime && x.charge == y.charge;
    }
  };

  struct DinosaurFeatureHash : std::unary_function<DinosaurFeature, std::size_t> {
    std::size_t operator()(const DinosaurFeature& x) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, x.precMz);
      boost::hash_combine(seed, x.rTime);
      boost::hash_combine(seed, x.charge);
      return seed;
    }
  };

 public:
  typedef FeatureList::iterator iterator;
  typedef FeatureList::const_iterator const_iterator;

  DinosaurFeatureList() : isInitialized_(false), features_() {};
  ~DinosaurFeatureList(){};

  /* delegation functions */
  inline iterator begin() { return features_.begin(); }
  inline iterator end() { return features_.end(); }
  inline const_iterator begin() const { return features_.begin(); }
  inline const_iterator end() const { return features_.end(); }
  inline std::vector<DinosaurFeature>& getFeatureList() { return features_; }
  inline const std::vector<DinosaurFeature>& getFeatureList() const { return features_; }
  inline size_t size() const { return features_.size(); }
  inline void clear() {
    std::vector<DinosaurFeature>().swap(features_);
    clearFeatureToIdxMap();
  }
  inline void push_back(const DinosaurFeature& ft) {
    features_.push_back(ft);
  }
  inline void push_back_with_ft_to_idx_map(const DinosaurFeature& ft) {
    features_.push_back(ft);
    featureToFeatureIdxMap_[ft] = ft.featureIdx;
  }
  inline const DinosaurFeature& at(size_t n) const { return features_.at(n); }
  inline DinosaurFeature& at(size_t n) { return features_.at(n); }

  inline bool isInitialized() const { return isInitialized_; }
  inline void setInitialized() { isInitialized_ = true; }

  inline void clearFeatureToIdxMap() { featureToFeatureIdxMap_.clear(); }
  inline int getFeatureIdx(const DinosaurFeature& ft) const {
    if (featureToFeatureIdxMap_.find(ft) != featureToFeatureIdxMap_.end()) {
      return featureToFeatureIdxMap_.find(ft)->second;
    } else {
      return -1;
    }
  }
  
  inline int getVectorIdx(size_t featureIdx) const {
    if (featureIdxToVectorIdxMap_.find(featureIdx) != featureIdxToVectorIdxMap_.end()) {
      return featureIdxToVectorIdxMap_.find(featureIdx)->second;
    } else {
      return -1;
    }
  }

  size_t loadFromFile(const std::string& ftFile, bool withIdxMap = false, 
      bool ignoreDuplicates = false, bool removePlaceholders = false) {
    std::vector<DinosaurFeature> addedFts;
    maracluster::BinaryInterface::read(ftFile, addedFts);
    features_.reserve(features_.size() + addedFts.size());
    for (DinosaurFeature& df : addedFts) {
      if ((!removePlaceholders || df.intensity > 0.0) && 
              (!ignoreDuplicates || getFeatureIdx(df) == -1)) { /* this mimicks the MBR feature adding */
        if (ignoreDuplicates) {
          df.featureIdx = size(); /* re-index before adding to the featurelist */
        }
        
        if (withIdxMap) {
          push_back_with_ft_to_idx_map(df);
        } else {
          push_back(df);
        }
      }
    }
    return addedFts.size();
  }

  void saveToFile(const std::string& ftFile, bool append) {
    maracluster::BinaryInterface::write(features_, ftFile, append);
  }
  
  void buildFeatureIdxToVectorIdxMap() {
    size_t featureIdx = 0u;
    for (const DinosaurFeature& df : features_) {
      featureIdxToVectorIdxMap_[df.featureIdx] = featureIdx++;
    }
  }

  inline void sortByPrecMz() { std::sort(features_.begin(), features_.end(), lessPrecMz); }
  inline void sortByFeatureIdx() { std::sort(features_.begin(), features_.end(), lessFeatureIdx); }

  inline static bool lessPrecMz(const DinosaurFeature& a, const DinosaurFeature& b) {
    return (a.precMz < b.precMz) || (a.precMz == b.precMz && a.rTime < b.rTime);
  }
  inline static bool lessFeatureIdx(const DinosaurFeature& a, const DinosaurFeature& b) { return (a.featureIdx < b.featureIdx); }
  inline std::vector<DinosaurFeature>::const_iterator getPrecMzIterator(float precMz) const {
    return std::lower_bound(begin(), end(), precMz, DinosaurFeatureLessPrecMzLower());
  }

 protected:
  FeatureList features_;
  bool isInitialized_;
  
  /* this map is used in FeatureAlignment.cpp to find features in the 
     DinosaurFeatureList based on a deserialized Percolator PSM id.
     See getFeatureIdx() */
  boost::unordered_map<DinosaurFeature, int, DinosaurFeatureHash, DinosaurFeatureEqual> featureToFeatureIdxMap_;
  boost::unordered_map<int, int> featureIdxToVectorIdxMap_;
};

} /* namespace quandenser */

#endif /* QUANDENSER_DINOSAURFEATURELIST_H_ */
