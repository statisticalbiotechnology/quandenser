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

#ifndef QUANDENSER_DINOSAURIO_H_
#define QUANDENSER_DINOSAURIO_H_

#include <vector>
#include <sstream>
#include <cctype>

#include "percolator/src/DataSet.h"

#include "Globals.h"
#include "MyException.h"
#include "DinosaurFeature.h"
#include "DinosaurFeatureList.h"

namespace quandenser {

class DinosaurIO {
 public:
  DinosaurIO(){};
  ~DinosaurIO(){};
  
  static void setJavaMemory(std::string mem);
  static void setJavaNumThreads(int numThreads) {
    javaNumThreads_ = numThreads;
  }
  static void setSeed(int seed) {
    seed_ = seed;
  }
  
  static void parseDinosaurFeatureFile(std::istream& dataStream, int fileIdx, DinosaurFeatureList& dinosaurFeatures);
  static int runDinosaurGlobal(const std::string& outputDir, const std::string& mzMLFile);
  static int runDinosaurTargeted(const std::string& outputDir, const std::string& mzMLFile, const std::string& targetFile);
  static DinosaurFeature parseDinosaurFeatureRow(std::string& line);
  
 protected:
  static int runDinosaur(const std::string& dinosaurCmd);
  static const std::string kAdvancedParamFile, kDinosaurJar;
  
  static std::string javaMemory_;
  static int javaNumThreads_;
  static int seed_;
  
    
};

} /* namespace quandenser */

#endif /* QUANDENSER_DINOSAURIO_H_ */
