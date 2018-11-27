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

#include "MaRaClusterAdapter.h"

namespace quandenser {

int MaRaClusterAdapter::createIndex() {
  if (spectrumBatchFileFN_.empty()) {
    std::cerr << "Error: no batch file specified with -b flag" << std::endl;
    return EXIT_FAILURE;
  }
  
  if (peakCountFN_.empty())
    peakCountFN_ = outputFolder_ + "/" + fnPrefix_ + ".peak_counts.dat";
  if (scanInfoFN_.empty())
    scanInfoFN_ = outputFolder_ + "/" + fnPrefix_ + ".scan_info.dat";
  if (datFNFile_.empty())
    datFNFile_ = outputFolder_ + "/" + fnPrefix_ + ".dat_file_list.txt";
  
  maracluster::SpectrumFileList fileList;
  fileList.initFromFile(spectrumBatchFileFN_);
  
  if (!Globals::fileExists(datFNFile_) || !Globals::fileExists(scanInfoFN_)) {    
    SpectrumFiles spectrumFiles(outputFolder_, chargeUncertainty_, 
                                featureLists_, spectrumToPrecursorMap_);
    spectrumFiles.splitByPrecursorMz(fileList, datFNFile_, peakCountFN_, 
        scanInfoFN_, precursorTolerance_, precursorToleranceDa_);
  } else {
    std::cerr << "Read dat-files from " << datFNFile_ << 
        " and scan numbers from " << scanInfoFN_ <<
        ". Remove these files to generate new dat-files." << std::endl;
  }
  
  return EXIT_SUCCESS;
}

int MaRaClusterAdapter::mergeSpectra() {
  if (spectrumOutFN_.empty())
    spectrumOutFN_ = outputFolder_ + "/" + fnPrefix_ + ".consensus.ms2";
  
  if (!maracluster::MSFileHandler::validMs2OutputFN(spectrumOutFN_)) {
    return EXIT_FAILURE;
  }
  
  if (maracluster::MSFileHandler::getOutputFormat(spectrumOutFN_) == "mgf") {
    maracluster::MSFileHandler::splitMassChargeStates_ = true;
  }
  
  ConsensusMerger msFileMerger(spectrumOutFN_, spectrumClusterToConsensusFeatures_);
  
  std::cerr << "Parsing cluster file" << std::endl;
  msFileMerger.parseClusterFileForMerge(clusterFileFN_, minConsensusClusterSize_);
  std::cerr << "Finished parsing cluster file" << std::endl;
  
  std::cerr << "Merging clusters" << std::endl;
  msFileMerger.mergeSpectra();
  std::cerr << "Finished merging clusters" << std::endl;
  
  return EXIT_SUCCESS;
}

} /* namespace quandenser */
