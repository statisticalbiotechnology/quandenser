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

#include "SpectrumFiles.h"

namespace quandenser {

void SpectrumFiles::getMassChargeCandidates(pwiz::msdata::SpectrumPtr s, 
    std::vector<maracluster::MassChargeCandidate>& mccs, maracluster::ScanId scanId) {
  //maracluster::SpectrumHandler::getMassChargeCandidates(s, mccs, chargeUncertainty_);
  //mccs.clear();
  if (!spectrumToPrecursorMap_.isInitialized(scanId)) {
    spectrumToPrecursorMap_.setInitialized(scanId);
    float precMz = maracluster::SpectrumHandler::getPrecMz(s);
    float rTime = maracluster::SpectrumHandler::getRetentionTime(s);  
    
    float isoWidthLower = s->precursors.at(0).isolationWindow.cvParam(pwiz::cv::MS_isolation_window_lower_offset).valueAs<float>();
    float isoWidthUpper = s->precursors.at(0).isolationWindow.cvParam(pwiz::cv::MS_isolation_window_upper_offset).valueAs<float>();
    
    // TODO: try if using a separate vector with only the precMzs is faster
    DinosaurFeatureList::const_iterator lowerBound = featureLists_.at(scanId.fileIdx).getPrecMzIterator(precMz - isoWidthLower);
    DinosaurFeatureList::const_iterator upperBound = featureLists_.at(scanId.fileIdx).getPrecMzIterator(precMz + isoWidthUpper);
    for (DinosaurFeatureList::const_iterator it = lowerBound; it != upperBound; ++it) {
      if (it->intensity > 0.0 && it->rtStart <= rTime && it->rtEnd >= rTime) {
        mccs.push_back(maracluster::MassChargeCandidate(it->charge, it->precMz, maracluster::SpectrumHandler::calcMass(it->precMz, it->charge)));
        DinosaurFeature featureCopy = *it;
        spectrumToPrecursorMap_.addFeature(scanId, featureCopy);
      }
    }
    
    if (Globals::VERB > 4 && mccs.size() == 0) {
      std::cerr << "could not find precursor: fileIdx = " << scanId.fileIdx << ", precmz = " << precMz << ", rtime = " << rTime << std::endl;
    }
  } else {
    for (DinosaurFeatureList::const_iterator it = spectrumToPrecursorMap_.getBegin(scanId); it != spectrumToPrecursorMap_.getEnd(scanId); ++it) {
      mccs.push_back(maracluster::MassChargeCandidate(it->charge, it->precMz, maracluster::SpectrumHandler::calcMass(it->precMz, it->charge)));
    }
  }
  
  if (mccs.size() > 1) {
    std::sort(mccs.begin(), mccs.end(), maracluster::MassChargeCandidate::lessChargeMass);
  }
}

} /* namespace quandenser */
