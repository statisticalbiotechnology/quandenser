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

#ifndef QUANDENSER_SPECTRUMFILES_H_
#define QUANDENSER_SPECTRUMFILES_H_

#include <vector>
#include <sstream>

#include "maracluster/src/SpectrumFiles.h"
#include "maracluster/src/SpectrumHandler.h"
#include "maracluster/src/MassChargeCandidate.h"
#include "maracluster/src/ScanId.h"

#include "Globals.h"
#include "DinosaurFeature.h"
#include "SpectrumToPrecursorMap.h"

namespace quandenser {

class SpectrumFiles : public maracluster::SpectrumFiles {
 public:
  SpectrumFiles(const std::string& precMzFileFolder, 
               const int chargeUncertainty, 
               const std::vector<DinosaurFeatureList>& featureLists,
               SpectrumToPrecursorMap& specToPrecMap) : 
    maracluster::SpectrumFiles(precMzFileFolder, chargeUncertainty), featureLists_(featureLists), spectrumToPrecursorMap_(specToPrecMap) {}
  
 protected:
  void getMassChargeCandidates(pwiz::msdata::SpectrumPtr s, 
    std::vector<maracluster::MassChargeCandidate>& mccs, maracluster::ScanId scanId);
  
  const std::vector<DinosaurFeatureList>& featureLists_;
  SpectrumToPrecursorMap& spectrumToPrecursorMap_;
    
};

} /* namespace quandenser */

#endif /* QUANDENSER_SPECTRUMFILES_H_ */
