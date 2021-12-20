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

#include "DinosaurIO.h"
#include <boost/lexical_cast.hpp>

namespace quandenser {

const std::string DinosaurIO::kAdvancedParamFileTargeted = "\"" + Globals::getJarPath() + std::string("advParams_dinosaur_targeted.txt") + "\"";
const std::string DinosaurIO::kAdvancedParamFile = "\"" + Globals::getJarPath() + std::string("advParams_dinosaur.txt") + "\"";
const std::string DinosaurIO::kDinosaurJar = "\"" + Globals::getJarPath() + std::string("Dinosaur-1.2.1.free.jar") + "\"";

std::string DinosaurIO::javaMemory_ = "24000M";
int DinosaurIO::javaNumThreads_ = 4;
int DinosaurIO::seed_ = 1;

void DinosaurIO::setJavaMemory(std::string mem) {
  char lastChar = std::tolower(*mem.rbegin());
  if ((lastChar == 'g' || lastChar == 'm') && boost::lexical_cast<int>(mem.substr(0, mem.size() - 1)) ) {
    javaMemory_ = mem;
  } else {
    std::ostringstream oss;
    oss << "Invalid memory string for Java's -Xmx argument: " << mem << ".\n"
        << "String should consist of a number followed by a \"G\" for gigabytes or \"M\" for megabytes.\n"
        << "Terminating.." << std::endl;
    throw MyException(oss.str());
  }
}

int DinosaurIO::runDinosaurGlobal(const std::string& outputDir, const std::string& mzMLFile) {
  return runDinosaur("--outDir=" + outputDir + " --advParams=" + kAdvancedParamFile + " " + mzMLFile);
}

int DinosaurIO::runDinosaurTargeted(const std::string& outputDir, const std::string& mzMLFile, const std::string& targetFile) {
  return runDinosaur("--outDir=" + outputDir + " --advParams=" + kAdvancedParamFileTargeted + " --mode=target --targets=" + targetFile + " " + mzMLFile);
}

int DinosaurIO::runDinosaur(const std::string& dinosaurCmd) {
  std::string cmd = "java -Xmx" + javaMemory_ + " -jar " + kDinosaurJar +
      " --force  --profiling=true --nReport=0 " +
      " --concurrency=" + boost::lexical_cast<std::string>(javaNumThreads_) +
      " --seed=" + boost::lexical_cast<std::string>(seed_) + " " + 
      dinosaurCmd;
  if (Globals::VERB > 2) {
    std::cerr << cmd << std::endl;
  }
  return std::system(cmd.c_str());
}

/*
  mz      mostAbundantMz  charge  rtStart rtApex  rtEnd   fwhm    nIsotopes       nScans  averagineCorr   mass    massCalib       intensityApex   intensitySum
  
  float precMz, rTime, intensity;
  int charge;
  int fileIdx, featureIdx;
  **/
void DinosaurIO::parseDinosaurFeatureFile(std::istream& dataStream, int fileIdx, DinosaurFeatureList& dinosaurFeatures) {
  std::string featureLine;
  
  getline(dataStream, featureLine); // skip header
  size_t lineNr = 2;
  std::vector<DinosaurFeature> ftVec;
  while (getline(dataStream, featureLine)) {
    if (lineNr % 1000000 == 0 && Globals::VERB > 1) {
      std::cerr << "Processing line " << lineNr << std::endl;
    }
    DinosaurFeature ft = parseDinosaurFeatureRow(featureLine);
    ftVec.push_back(ft);
    
    ++lineNr;
  }
  
  std::sort(ftVec.begin(), ftVec.end(), DinosaurFeatureList::lessPrecMz);
  
  size_t ftIdx = 0;
  for (std::vector<DinosaurFeature>::iterator it = ftVec.begin(); it != ftVec.end(); ++it, ++ftIdx) {
    it->featureIdx = ftIdx;
    it->fileIdx = fileIdx;
    
    dinosaurFeatures.push_back(*it);
  }
}

DinosaurFeature DinosaurIO::parseDinosaurFeatureRow(std::string& line) {
  TabReader reader(line);
  
  DinosaurFeature ft;
  
  ft.precMz = static_cast<float>(reader.readDouble()); // mz
  reader.readDouble(); // mostAbundantMz
  ft.charge = reader.readInt(); // charge
  ft.rtStart = static_cast<float>(reader.readDouble()); // rtStart
  ft.rTime = static_cast<float>(reader.readDouble()); // rtApex
  ft.rtEnd = static_cast<float>(reader.readDouble()); // rtEnd
  reader.readDouble(); // fwhm
  reader.readInt(); // nIsotopes
  reader.readInt(); // nScans
  reader.readDouble(); // averagineCorr
  reader.readDouble(); // mass
  reader.readDouble(); // massCalib
  reader.readDouble(); // intensityApex
  ft.intensity = static_cast<float>(reader.readDouble()); // intensitySum
  
  return ft;
}

} /* namespace quandenser */
