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

#ifndef QUANDENSER_PERCOLATORADAPTER_H_
#define QUANDENSER_PERCOLATORADAPTER_H_

#include <vector>
#include <map>
#include <sstream>

#include "percolator/src/Caller.h"

#include "Globals.h"
#include "DinosaurFeature.h"

namespace quandenser {

class PercolatorAdapter : public Caller {
 public:
  PercolatorAdapter();
  
  bool parseOptions(const std::vector<std::string>& percolatorArgs) {
    std::vector<const char*> percolatorArgv;
    for (std::vector<std::string>::const_iterator it = percolatorArgs.begin();
       it != percolatorArgs.end(); ++it) {
      percolatorArgv.push_back(it->c_str());
    }
    return Caller::parseOptions(percolatorArgs.size(), (char**)&percolatorArgv.front());
  }
  
  void init(const std::vector<std::pair<std::string, double> >& percolatorFeatures);
  void addPsm(const std::string& specId, int label, int scannr, 
    const std::string& peptide, const std::vector<double>& features);
  
  void process();
 protected:
  SetHandler setHandler_;
  
};

} /* namespace quandenser */

#endif /* QUANDENSER_PERCOLATORADAPTER_H_ */
