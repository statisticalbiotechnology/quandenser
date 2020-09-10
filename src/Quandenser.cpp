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

#include "Quandenser.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = boost::filesystem;

namespace quandenser {

Quandenser::Quandenser() : call_(""), fnPrefix_("Quandenser"), seed_(1u),
    numThreads_(4u), maxMissingValues_(-1), intensityScoreThreshold_(0.5),
    spectrumBatchFileFN_(""), outputFolder_("Quandenser_output"),
    outputSpectrumFile_(""), maraclusterPpmTol_(20.0f),
    alignPpmTol_(20.0f), alignRTimeStdevTol_(10.0f),
    decoyOffset_(5.0 * 1.000508),
    linkPEPThreshold_(0.25), linkPEPMbrSearchThreshold_(0.05),
    maxFeatureCandidates_(2),  useTempFiles_(false),
    partial1Dinosaur_(-1), partial2MaRaCluster_(""), 
    partial3MatchRound_(-1), partial4Consensus_("") {}

std::string Quandenser::greeter() {
  std::ostringstream oss;
  oss << "Quandenser version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Copyright (c) 2018-19 Matthew The. All rights reserved.\n"
      << "Written by Matthew The (matthew.the@scilifelab.se) in the\n"
      << "School of Biotechnology at the Royal Institute of Technology in Stockholm.\n";
  return oss.str();
}

std::string Quandenser::extendedGreeter(time_t& startTime) {
  std::ostringstream oss;
  char* host = getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << endl << call_ << endl;
  oss << "Started " << ctime(&startTime) << endl;
  oss.seekp(-1, ios_base::cur);
  if (host) oss << " on " << host << endl;
  return oss.str();
}

bool Quandenser::parseOptions(int argc, char **argv) {
  std::ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << endl;
  call_ = callStream.str();
  call_ = call_.substr(0,call_.length()-1); // trim ending carriage return

  maraclusterArgs_.push_back("maracluster");
  maraclusterArgs_.push_back("mode_placeholder");
  maraclusterArgs_.push_back("--splitMassChargeStates");

  percolatorArgs_.push_back("percolator");
  percolatorArgs_.push_back("--only-psms");
  percolatorArgs_.push_back("--post-processing-tdc");
  /* prevents percolator's Caller::parseOptions from throwing an error because
     no input file was specified */
  percolatorArgs_.push_back("input_file_placeholder");

  std::ostringstream intro;
  intro << greeter() << "\nUsage:\n" << endl;
  intro << "  quandenser -b <batch_file_list> -f <output_directory>" << std::endl;

  // init
  // available 1-letter abbreviations: d,e,g,i,j,k,q,s,u,v,w,x,y,z
  // A,C,D,E,F,G,H,J,K,L,N,O,Q,R,S,T,U,V,W,X,Y,Z
  CommandLineParser cmd(intro.str());
  cmd.defineOption("b",
      "batch",
      "File with spectrum files to be processed in batch, one per line",
      "filename");
  cmd.defineOption("f",
      "output-folder",
      "Writable folder for output files (default: 'Quandenser_output')",
      "path");
  cmd.defineOption("s",
      "seed",
      "Seed for pseudo random number generators of Dinosaur and Percolator (default: 1)",
      "int");
  cmd.defineOption("a",
      "prefix",
      "Output files will be prefixed as e.g. <prefix>.clusters_p10.tsv (default: 'Quandenser')",
      "name");
  cmd.defineOption("m",
      "max-missing",
      "Maximum number of missing values in the quantification row (default: numFiles / 4)",
      "name");
  cmd.defineOption("v",
      "verbatim",
      "Set the verbatim level (lowest: 0, highest: 5, default: 3).",
      "int");
  cmd.defineOption("N",
      "num-threads",
      "Number of threads used for Quandenser, this includes the number of threads for MaRaCluster, Percolator and Dinosaur (the latter can be overriden by --dinosaur-threads) (default: 4).",
      "int");
  cmd.defineOption("o",
      "spec-out",
      "Path to output file for the consensus spectra, the infix \".part<x>\" will be added before the extension. All output file formats supported by ProteoWizard are valid (default: <output-folder>/consensus_spectra/<prefix>.consensus.ms2)",
      "filename");
  cmd.defineOption("t",
      "percolator-test-fdr",
      "False discovery rate threshold for evaluating best cross validation result and reported end result (default: 0.02).",
      "value");
  cmd.defineOption("F",
      "percolator-train-fdr",
      "False discovery rate threshold to define positive examples in training. Set to testFDR if 0 (default: 0.02).",
      "value");
  cmd.defineOption("M",
      "dinosaur-memory",
      "Memory available for Dinosaur. String should consist of a number followed by a \"G\" for gigabytes or \"M\" for megabytes (default: 24000M).",
      "value");
  cmd.defineOption("n",
      "dinosaur-threads",
      "Maximum number of threads available for Dinosaur, will override --num-threads if specified (default: 4).",
      "int");
  cmd.defineOption("p",
      "maracluster-pval-threshold",
      "Set log(p-value) threshold for MaRaCluster MS/MS clustering (default: -10.0).",
      "double");
  cmd.defineOption("z",
      "maracluster-mz-tol",
      "Mass tolerance for clustering fragment spectra with MaRaCluster in ppm (default: 20.0).",
      "double");
  cmd.defineOption("P",
      "align-mz-tol",
      "Mass tolerance for matching features between runs in ppm (default: 20.0).",
      "double");
  cmd.defineOption("r",
      "align-rtime-tol",
      "Retention time tolerance for matching feature between runs in number of standard deviations of the retention time alignment (default: 10.0).",
      "double");
  cmd.defineOption("I",
      "intensity-score-cut",
      "Intensity score ratio cut-off relative to the highest intensity score per cluster (lowest: 0.0, highest: 1.0, default: 0.5).",
      "double");
  cmd.defineOption("l",
      "ft-link-cut",
      "Posterior error probability cut-off for matching features between runs (lowest: 0.0, highest: 1.0, default: 0.25).",
      "double");
  cmd.defineOption("c",
      "ft-link-candidates",
      "Number of candidates to consider for matching features between runs (default: 2).",
      "int");
  cmd.defineOption("B",
      "target-search-threshold",
      "Minimum posterior error probability for re-searching for a feature with a targeted search. Setting this to 1.0 will cause the targeted search to be skipped (lowest: 0.0, highest: 1.0, default: 0.05).",
      "double");
  cmd.defineOption("W",
      "partial-1-dinosaur",
      "Partial run step 1. Runs the initial Dinosaur feature finding runs. Specify either the file index as specified in the batch input file (starting with file index 1), or specify 0 to process all files in one go.",
      "int");
  cmd.defineOption("X",
      "partial-2-maracluster",
      "Partial run step 2. Runs the first MaRaCluster run and the retention time alignments.",
      "mode");
  cmd.defineOption("Y",
      "partial-3-match-round",
      "Partial run step 3. Runs a feature matching round from the minimum spanning tree of run alignments. Specify either the round number as specified in the feature alignment queue (starting at round 1), or specify 0 to process all rounds in one go.",
      "int");
  cmd.defineOption("Z",
      "partial-4-consensus",
      "Partial run step 4. Runs second MaRaCluster run, creates and write feature groups and consensus spectra.",
      "mode");
  cmd.defineOption(Option::EXPERIMENTAL_FEATURE,
      "decoy-offset",
      "Decoy m/z offset (lowest: -1000.0, highest: 1000.0, default: 5.0 * 1.000508).",
      "double");
  cmd.defineOption(Option::EXPERIMENTAL_FEATURE,
      "use-tmp-files",
      "Write temporary files to disk at several points in the process to reduce memory usage (default: false).",
      "",
      TRUE_IF_SET);
  /*
    maxFeatureCandidates_(2)
    */
  // parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);

  if (cmd.optionSet("batch")) {
    spectrumBatchFileFN_ = cmd.options["batch"];

    maraclusterArgs_.push_back("--batch");
    maraclusterArgs_.push_back(cmd.options["batch"]);
  }

  if (cmd.optionSet("output-folder")) {
    outputFolder_ = cmd.options["output-folder"];
  }

  if (cmd.optionSet("seed")) {
    int seed = cmd.getInt("seed", 1, 20000); // [1,20000] are Percolator's limits for the seed
    percolatorArgs_.push_back("--seed");
    percolatorArgs_.push_back(cmd.options["seed"]);
    DinosaurIO::setSeed(seed);
  }

  if (cmd.optionSet("spec-out")) {
    outputSpectrumFile_ = cmd.options["spec-out"];
  }

  if (cmd.optionSet("prefix")) {
    fnPrefix_ = cmd.options["prefix"];

    maraclusterArgs_.push_back("--prefix");
    maraclusterArgs_.push_back(cmd.options["prefix"]);
  }

  if (cmd.optionSet("max-missing")) {
    maxMissingValues_ = cmd.getInt("max-missing", 0, 1000);;
  }

  if (cmd.optionSet("verbatim")) {
    int v = cmd.getInt("verbatim", 0, 5);
    Globals::VERB = v;

    maraclusterArgs_.push_back("--verbatim");
    maraclusterArgs_.push_back(cmd.options["verbatim"]);

    percolatorArgs_.push_back("--verbatim");
    percolatorArgs_.push_back(cmd.options["verbatim"]);
  }

  if (cmd.optionSet("num-threads")) {
    numThreads_ = cmd.getInt("num-threads", 1, 1000);
    percolatorArgs_.push_back("--num-threads");
    percolatorArgs_.push_back(cmd.options["num-threads"]);
    DinosaurIO::setJavaNumThreads(cmd.getInt("num-threads", 1, 1000));
  }

  if (cmd.optionSet("dinosaur-memory")) {
    DinosaurIO::setJavaMemory(cmd.options["dinosaur-memory"]);
  }

  if (cmd.optionSet("dinosaur-threads")) {
    DinosaurIO::setJavaNumThreads(cmd.getInt("dinosaur-threads", 1, 1000));
  }

  if (cmd.optionSet("maracluster-pval-threshold")) {
    /* just to test if value is valid, MaRaCluster will parse the string */
    float pvalThresholdTest = cmd.getDouble("maracluster-pval-threshold", -1000.0, 0.0);
    maraclusterArgs_.push_back("--clusterThresholds");
    maraclusterArgs_.push_back(cmd.options["maracluster-pval-threshold"]);
    maraclusterArgs_.push_back("--pvalThreshold");
    maraclusterArgs_.push_back(cmd.options["maracluster-pval-threshold"]);
  } else {
    maraclusterArgs_.push_back("--clusterThresholds");
    maraclusterArgs_.push_back("-10.0");
    maraclusterArgs_.push_back("--pvalThreshold");
    maraclusterArgs_.push_back("-10.0");
  }

  percolatorArgs_.push_back("--trainFDR");
  if (cmd.optionSet("percolator-train-fdr")) {
    /* just to test if value is valid, MaRaCluster will parse the string */
    float trainFdrTest = cmd.getDouble("percolator-train-fdr", 0.0, 1.0);
    percolatorArgs_.push_back(cmd.options["percolator-train-fdr"]);
  } else {
    percolatorArgs_.push_back("0.02");
  }

  percolatorArgs_.push_back("--testFDR");
  if (cmd.optionSet("percolator-test-fdr")) {
    float testFdrTest = cmd.getDouble("percolator-test-fdr", 0.0, 1.0);
    percolatorArgs_.push_back(cmd.options["percolator-train-fdr"]);
  } else {
    percolatorArgs_.push_back("0.02");
  }

  if (cmd.optionSet("maracluster-mz-tol")) {
    MaRaClusterIO::setPrecursorTolerance(cmd.getDouble("maracluster-mz-tol", 0.0, 100000.0));
  }

  if (cmd.optionSet("align-mz-tol")) {
    alignPpmTol_ = cmd.getDouble("align-mz-tol", 0.0, 100000.0);
  }

  if (cmd.optionSet("align-rtime-tol")) {
    alignRTimeStdevTol_ = cmd.getDouble("align-rtime-tol", 0.0, 1000.0);
  }

  if (cmd.optionSet("intensity-score-cut")) {
    intensityScoreThreshold_ = cmd.getDouble("intensity-score-cut", 0.0, 1.0);
  }

  if (cmd.optionSet("ft-link-cut")) {
    linkPEPThreshold_ = cmd.getDouble("ft-link-cut", 0.0, 1.0);
  }

  if (cmd.optionSet("ft-link-candidates")) {
    linkPEPThreshold_ = cmd.getInt("ft-link-candidates", 1, 100);
  }

  if (cmd.optionSet("target-search-threshold")) {
    maxFeatureCandidates_ = cmd.getDouble("target-search-threshold", 0.0, 1.0);
  }

  if (cmd.optionSet("partial-1-dinosaur")) {
    partial1Dinosaur_ = cmd.getInt("partial-1-dinosaur", 0, 99999);
  }

  if (cmd.optionSet("partial-2-maracluster")) {
    partial2MaRaCluster_ = cmd.options["partial-2-maracluster"];
  }

  if (cmd.optionSet("partial-3-match-round")) {
    partial3MatchRound_ = cmd.getInt("partial-3-match-round", 0, 99999);
  }

  if (cmd.optionSet("partial-4-consensus")) {
    partial4Consensus_ = cmd.options["partial-4-consensus"];
  }

  if (cmd.optionSet("decoy-offset")) {
    decoyOffset_ = cmd.getDouble("decoy-offset", -1000.0, 1000.0);
  }

  if (cmd.optionSet("use-tmp-files")) {
    useTempFiles_ = true;
  }

  // if there are arguments left...
  if (cmd.arguments.size() > 0) {
    std::cerr << "Error: too many arguments." << std::endl;
    std::cerr << "Invoke with -h option for help." << std::endl;
    return 0; // ...error
  } else {
    if (!cmd.optionSet("batch")) {
      std::cerr << "Error: one of the inputs is missing." << std::endl;
      std::cerr << "Invoke with -h option for help." << std::endl;
      return 0;
    }
  }

  return true;
}

int Quandenser::runMaRaCluster(const std::string& maRaClusterSubFolder,
    const maracluster::SpectrumFileList& fileList,
    std::vector<DinosaurFeatureList>& allFeatures,
    std::string& clusterFilePath,
    const std::string& spectrumToPrecursorFile,
    const std::string& tmpFilePrefixAlign,
    const std::string& mode) {
  std::vector<std::string> maraclusterArgs = maraclusterArgs_;
  fs::path maraclusterFolder(outputFolder_);
  maraclusterFolder /= maRaClusterSubFolder;
  if (mode == "index") {
    maraclusterArgs[1] = "index";
  } else if (mode.rfind("pvalue", 0) == 0) {
    maraclusterArgs[1] = "pvalue";
    std::string datFile = mode.substr(7);
    maraclusterArgs.push_back("--specIn");
    maraclusterArgs.push_back(maraclusterFolder.string() + "/" + datFile);
    maraclusterArgs.push_back("--peakCountsFN");
    maraclusterArgs.push_back(maraclusterFolder.string() + "/MaRaCluster.peak_counts.dat");
    maraclusterArgs.push_back("--prefix");
    maraclusterArgs.push_back(datFile);
    maraclusterArgs.push_back("--clusteringTree");
    maraclusterArgs.push_back(maraclusterFolder.string() + "/" + datFile + ".pvalue_tree.tsv");
    maraclusterArgs.push_back("--pvalOut");
    maraclusterArgs.push_back(maraclusterFolder.string() + "/" + datFile + ".pvalues.dat");    
  } else if (mode.rfind("overlap", 0) == 0) {
    maraclusterArgs[1] = "overlap";
    std::string overlapIdx = mode.substr(8);
    maraclusterArgs.push_back("--datFNfile");
    maraclusterArgs.push_back(maraclusterFolder.string() + "/MaRaCluster.dat_file_list.txt");
    maraclusterArgs.push_back("--overlapBatchIdx");
    maraclusterArgs.push_back(overlapIdx);
  } else {
    maraclusterArgs[1] = "batch";
  }
  maraclusterArgs.push_back("--output-folder");
  maraclusterArgs.push_back(maraclusterFolder.string());
  maraclusterArgs.push_back("--scanInfoFN");
  maraclusterArgs.push_back(maraclusterFolder.string() + "/MaRaCluster.scan_info.dat");

  std::string spectrumOutputFile = ""; /* Not needed here */
  SpectrumToPrecursorMap spectrumToPrecursorMap(fileList.size());
  MaRaClusterAdapter maraclusterAdapter(allFeatures, spectrumToPrecursorMap, spectrumOutputFile);
  maraclusterAdapter.parseOptions(maraclusterArgs);
  
  clusterFilePath = maraclusterAdapter.getClusterFileName();
  
  bool createIndex = (mode.empty() || mode == "index" || mode == "batch");
  if (!fs::exists(spectrumToPrecursorFile) || !mode.empty()) {
    if (useTempFiles_ && createIndex) {
      loadAllFeatures(tmpFilePrefixAlign, allFeatures);
    }
    
    if (!(mode == "cluster" && fs::exists(clusterFilePath))) {
      int rc = maraclusterAdapter.run();
      if (rc != EXIT_SUCCESS) {
        std::cerr << "MaRaCluster failed with exit code " << rc << ". Terminating.." << std::endl;
        return rc;
      }
    } else if (Globals::VERB > 1) {
      std::cerr << "Found cluster result file at " << clusterFilePath << ". Remove this file to redo clustering." << std::endl;
    }
      

    if (useTempFiles_ && createIndex) {
      unloadAllFeatures(allFeatures);
    }
    
    if (createIndex) {
      spectrumToPrecursorMap.serialize(spectrumToPrecursorFile);
      if (Globals::VERB > 1) {
        std::cerr << "Serialized spectrum to precursor map" << std::endl;
      }
    }
  }

  return EXIT_SUCCESS;
}

int Quandenser::saveAlignmentQueue(
    std::vector<std::pair<int, FilePair> >& featureAlignmentQueue,
    const std::string& featureAlignQueueFile,
    maracluster::SpectrumFileList& fileList) {
  std::ofstream outfile(featureAlignQueueFile.c_str(), ios::out);
  if (outfile.is_open()) {
    std::vector<std::pair<int, FilePair> >::const_iterator filePairIt;
    for (filePairIt = featureAlignmentQueue.begin();
          filePairIt != featureAlignmentQueue.end(); ++filePairIt) {
      FilePair filePair = filePairIt->second;
      int round = filePairIt->first + 1;
      fs::path file1(fileList.getFilePath(filePair.fileIdx1));
      fs::path file2(fileList.getFilePath(filePair.fileIdx2));
      // Serialized output stream
      outfile << round << '\t' << file1.string() << '\t' << file2.string()
              << '\t' << filePair.fileIdx1 << '\t' << filePair.fileIdx2 << std::endl;
    }
    outfile.close();
    return EXIT_SUCCESS;
  } else {
    std::cerr << "ERROR: Could not open " << featureAlignQueueFile
        << " for writing the feature alignment queue." << std::endl;
    return EXIT_FAILURE;
  }
}

std::string Quandenser::getFeatureFN(
    const fs::path& basePath, size_t fileIdx) {
  return FeatureAlignment::getFeatureFN(basePath.string(), fileIdx);
}

std::string Quandenser::getAddedFeaturesFN(
    const fs::path& basePath,
    size_t fileIdx1, size_t fileIdx2) {
  return FeatureAlignment::getAddedFeaturesFN(basePath.string(), fileIdx1, fileIdx2);
}

void Quandenser::loadAllFeatures(const std::string& tmpFilePrefixAlign,
    std::vector<DinosaurFeatureList>& allFeatures) {
  for (size_t i = 0; i < allFeatures.size(); ++i) {
    std::string fileName = getFeatureFN(fs::path(tmpFilePrefixAlign) / "features", i);
    bool withIdxMap = false;
    bool ignoreDuplicates = false;
    bool removePlaceholders = true;
    allFeatures.at(i).loadFromFile(fileName, withIdxMap, ignoreDuplicates, removePlaceholders);
    if (Globals::VERB > 2) {
      std::cerr << "Read in " << allFeatures.at(i).size() << " features from " << fileName << std::endl;
    }
  }
}

void Quandenser::unloadAllFeatures(std::vector<DinosaurFeatureList>& allFeatures) {
  for (size_t i = 0; i < allFeatures.size(); ++i) {
    allFeatures.at(i).clear();
  }
}

void Quandenser::updateMatchesTmpFile(const fs::path& matchesTmpFileBase, 
    const FilePair& filePair, 
    const DinosaurFeatureList& currentFeatures, 
    const DinosaurFeatureList& addedFeatures,
    const size_t originalNumFeaturesFromFile) {          
  std::string matchesFileName = getAddedFeaturesFN(
      matchesTmpFileBase, filePair.fileIdx1, filePair.fileIdx2);
  std::string matchesFileNameBackup = matchesFileName + ".bak";
  if (!fs::exists(matchesFileNameBackup)) {
    rename(matchesFileName.c_str(), matchesFileNameBackup.c_str());
    std::map<int, FeatureIdxMatch> featureMatches;
    FeatureAlignment::loadFromFile(matchesFileNameBackup, featureMatches);
    
    std::map<int, FeatureIdxMatch>::iterator ftIdxMatchIt;
    for (ftIdxMatchIt = featureMatches.begin(); 
         ftIdxMatchIt != featureMatches.end();
         ++ftIdxMatchIt) {
      int targetFeatureIdx = ftIdxMatchIt->second.targetFeatureIdx;
      if (targetFeatureIdx >= originalNumFeaturesFromFile) {
        ftIdxMatchIt->second.targetFeatureIdx = 
          currentFeatures.getFeatureIdx(addedFeatures.at(targetFeatureIdx - originalNumFeaturesFromFile));
      }
    }
    
    bool append = false;
    FeatureAlignment::saveToFile(matchesFileName, featureMatches, append);
  }
}

int Quandenser::createDirectory(fs::path dirPath) {
  boost::system::error_code returnedError;
  
  fs::create_directories(dirPath, returnedError);
  if (!fs::exists(dirPath)) {
    std::cerr << "Error: could not create output directory at "
        << dirPath.string() << std::endl;
    return EXIT_FAILURE;
  } else {
    return EXIT_SUCCESS;
  }
}

// adapted from https://github.com/crux-toolkit/crux-toolkit/blob/master/src/util/crux-utils.cpp
bool Quandenser::parseUrl(std::string url, std::string* host, std::string* path) {
  if (!host || !path) {
    return false;
  }
  // find protocol
  size_t protocolSuffix = url.find("://");
  if (protocolSuffix != std::string::npos) {
    url = url.substr(protocolSuffix + 3);
  }
  size_t pathBegin = url.find('/');
  if (pathBegin == std::string::npos) {
    *host = url;
    *path = "/";
  } else {
    *host = url.substr(0, pathBegin);
    *path = url.substr(pathBegin);
  }
  if (host->empty()) {
    *host = *path = "";
    return false;
  }
  return true;
}

void Quandenser::httpRequest(const std::string& url, const std::string& data) {
  // Parse URL into host and path components
  std::string host, path;
  if (!parseUrl(url, &host, &path)) {
    if (Globals::VERB > 2) {
      std::cerr << "Warning: Failed parsing URL " << url << std::endl;
    }
    return;
  }

  using namespace boost::asio;

  // Establish TCP connection to host on port 80
  io_service service;
  ip::tcp::resolver resolver(service);
  ip::tcp::resolver::iterator endpoint = resolver.resolve(ip::tcp::resolver::query(host, "80"));
  ip::tcp::socket sock(service);
  connect(sock, endpoint);
  
  std::size_t seed = 0;
  boost::hash_combine(seed, ip::host_name());
  boost::hash_combine(seed, sock.local_endpoint().address().to_string());
  std::stringstream stream;
  stream << std::hex << seed;
  
  std::string placeholder = "CID_PLACEHOLDER";
  std::string cid = stream.str();
  
  std::string newData(data);
  
  if (Globals::VERB > 3) {
    std::cerr << "Analytics data string: " << newData << std::endl;
  }
  
  newData.replace(newData.find(placeholder), placeholder.length(), cid);
  
  // Determine method (GET if no data; otherwise POST)
  std::string method = newData.empty() ? "GET" : "POST";
  std::ostringstream lengthString;
  lengthString << newData.length();
  
  std::string contentLengthHeader = newData.empty()
    ? ""
    : "Content-Length: " + lengthString.str() + "\r\n";
  // Send the HTTP request
  std::string request =
    method + " " + path + " HTTP/1.1\r\n"
    "Host: " + host + "\r\n" +
    contentLengthHeader +
    "Connection: close\r\n"
    "\r\n" + newData;
  sock.send(buffer(request));
}

void Quandenser::postToAnalytics(const std::string& appName) {
  // Post data to Google Analytics
  // For more information, see: https://developers.google.com/analytics/devguides/collection/protocol/v1/devguide
  try {
    std::stringstream paramBuilder;
    paramBuilder << "v=1"                // Protocol verison
                 << "&tid=UA-165948942-2" // Tracking ID
                 << "&cid=CID_PLACEHOLDER" // Unique device ID
                 << "&t=event"           // Hit type
                 << "&ec=quandenser"     // Event category
                 << "&ea=" << appName    // Event action
                 << "&el="               // Event label
#ifdef _MSC_VER
                      "win"
#elif __APPLE__
                      "mac"
#else
                      "linux"
#endif
                   << '-' << VERSION;
    httpRequest(
      "http://www.google-analytics.com/collect",
      paramBuilder.str());
  } catch (...) {
  }
}

int Quandenser::run() {
  time_t startTime;
  clock_t startClock;
  time(&startTime);
  startClock = clock();

  if (Globals::VERB > 0) {
    std::cerr << extendedGreeter(startTime);
  }

  std::string appName = "quandenser";
  postToAnalytics(appName);

#ifdef _OPENMP
  omp_set_num_threads(std::min((unsigned int)omp_get_max_threads(), numThreads_));
#endif

  boost::system::error_code returnedError;

  fs::path rootPath(outputFolder_);
  int rc = createDirectory(rootPath);
  if (rc != EXIT_SUCCESS) return rc;
    
  bool isPartialRun = (partial1Dinosaur_ >= 0 || !partial2MaRaCluster_.empty()
                  || partial3MatchRound_ >= 0 || !partial4Consensus_.empty());

  /*****************************************************************************
    Step 0: Read in the input text file with a list of mzML files
   ****************************************************************************/

  maracluster::SpectrumFileList fileList;
  fileList.initFromFile(spectrumBatchFileFN_);

  if (fileList.size() < 2u && partial1Dinosaur_ <= 0) {
    std::cerr << "Error: less than 2 spectrum files were specified, "
        << "Quandenser needs at least two files to perform a "
        << "meaningful alignment." << std::endl;
    return EXIT_FAILURE;
  }
  
  if (maxMissingValues_ < 0) {
    maxMissingValues_ = fileList.size() / 4;
  }

  /*****************************************************************************
    Step 1: Detect MS1 features with Dinosaur
   ****************************************************************************/

  fs::path dinosaurFolder(rootPath);
  dinosaurFolder /= "dinosaur";
  rc = createDirectory(dinosaurFolder);
  if (rc != EXIT_SUCCESS) return rc;
  
  /* binary feature file base path */
  fs::path featureOutFile(dinosaurFolder);
  featureOutFile /= "features";
  
  std::string tmpFilePrefixAlign = "";
  if (useTempFiles_) {
    fs::path tmpFileFolder(rootPath);
    tmpFileFolder /= "tmp";
    tmpFileFolder /= "matchFeatures";
    rc = createDirectory(tmpFileFolder);
    if (rc != EXIT_SUCCESS) return rc;
    
    tmpFilePrefixAlign = tmpFileFolder.string();
  }

  std::vector<std::string> files = fileList.getFilePaths();
  if (partial1Dinosaur_ > 0) {
    files = std::vector<std::string>(1, files.at(partial1Dinosaur_ - 1));
  }
  std::vector<fs::path> dinosaurFeatureFiles;
  std::vector<std::string>::const_iterator it = files.begin();
  for ( ; it != files.end(); ++it) {
    if (maracluster::MSFileHandler::getOutputFormat(*it) != "mzml" && Globals::VERB > 0) {
      std::cerr << "Warning: file extension is not .mzML, ignoring file: " << *it << std::endl;
    }
    fs::path mzMLFile(*it);
    fs::path dinosaurFeatureFile(dinosaurFolder.string());

    /* TODO: what if we use files with the same filename but in different folders? */
    dinosaurFeatureFile /= mzMLFile.stem().string() + ".features.tsv";
    dinosaurFeatureFiles.push_back(dinosaurFeatureFile);
    if (!isPartialRun || partial1Dinosaur_ >= 0) {
      if (!fs::exists(dinosaurFeatureFile)) {
        if (Globals::VERB > 1) {
          std::cerr << "Processing " << mzMLFile.filename() << " with Dinosaur." << std::endl;
        }
        rc = DinosaurIO::runDinosaurGlobal(dinosaurFolder.string(), mzMLFile.string());
        if (rc != EXIT_SUCCESS) {
          std::cerr << "Dinosaur failed with exit code " << rc << ". Terminating.." << std::endl;
          return EXIT_FAILURE;
        }
      } else if (Globals::VERB > 1) {
        std::cerr << "Already processed " << mzMLFile.filename() << " with Dinosaur." << std::endl;
      }
    }
  }

  std::vector<DinosaurFeatureList> allFeatures;
  if (!isPartialRun || partial1Dinosaur_ >= 0) {
    size_t fileIdx = 0u;
    for (std::vector<fs::path>::const_iterator it = dinosaurFeatureFiles.begin(); it != dinosaurFeatureFiles.end(); ++it, ++fileIdx) {
      size_t featureFileIdx = (partial1Dinosaur_ <= 0) ? fileIdx : partial1Dinosaur_ - 1;
      
      std::ifstream fileStream(it->string().c_str(), ios::in);
      DinosaurFeatureList features;
      DinosaurIO::parseDinosaurFeatureFile(fileStream, featureFileIdx, features);

      features.sortByPrecMz();

      allFeatures.push_back(features);
      if (Globals::VERB > 2) {
        std::cerr << "Read in " << features.size() << " features from " << it->filename() << std::endl;
      }
      
      if (partial1Dinosaur_ >= 0) {
        if (Globals::VERB > 2) {
          std::cerr << "Partial run: Saving dinosaur features in binary format for run " << featureFileIdx << std::endl;
        }
        // these files will not be overwritten by featureAlignment.matchFeatures
        std::string featureListFile = getFeatureFN(featureOutFile, featureFileIdx);
        bool append = false;
        features.saveToFile(featureListFile, append);
      }
      
      if (useTempFiles_) {
        std::string fileName = getFeatureFN(fs::path(tmpFilePrefixAlign) / "features", featureFileIdx);
        // remove the temp file, as it could have been overwritten by featureAlignment.matchFeatures
        remove(fileName.c_str()); 
        bool append = false;
        features.saveToFile(fileName, append);
        features.clear();
      }
    }
  }
  
  if (partial1Dinosaur_ >= 0) {
    if (Globals::VERB > 2) {
      std::cerr << "Partial run for Dinosaur feature detection completed" << std::endl;
    }
    return EXIT_SUCCESS;
  }
  
  if (isPartialRun) {
    if (Globals::VERB > 2) {
      std::cerr << "Partial run: Loading dinosaur features" << std::endl;
    }
    size_t fileIdx = 0u;
    std::vector<fs::path>::const_iterator it = dinosaurFeatureFiles.begin();
    for (; it != dinosaurFeatureFiles.end(); ++it, ++fileIdx) {
      DinosaurFeatureList ftList;
      allFeatures.push_back(ftList);
      fs::path ftFilePath(fileList.getFilePath(fileIdx));
      std::string fn(ftFilePath.filename().string());
      
      bool createMaRaClusterIndex = (!partial2MaRaCluster_.empty() &&
                                      (partial2MaRaCluster_ == "index" 
                                      || partial2MaRaCluster_ == "batch"))
                                     || (!partial4Consensus_.empty() &&
                                      (partial4Consensus_ == "index" 
                                      || partial4Consensus_ == "batch"
                                      || partial4Consensus_ == "cluster"));
      
      // Check if the pair assigned exists in the pair folder
      bool processRunInAlignment = (partial3MatchRound_ == 0
                                || (partial3MatchRound_ > 0 &&
                                     (fs::exists("pair/file1/" + fn)
                                      || fs::exists("pair/file2/" + fn)))); 
      
      /**
        If useTempFiles is set and partial runs are activated, we load the 
        features from temporary files. Otherwise, we need to get the features
        only in "batch" and "index" modes, as well as in the "cluster" mode
        for partial4Consensus_.
       **/
      if ((createMaRaClusterIndex && !useTempFiles_) || processRunInAlignment) {
        std::string featureListFile = getFeatureFN(featureOutFile, fileIdx);
        bool withIdxMap = processRunInAlignment;
        size_t ftsAdded = allFeatures.at(fileIdx).loadFromFile(featureListFile, withIdxMap);
        if (Globals::VERB > 2) {
          std::cerr << "Read in " << allFeatures.at(fileIdx).size() << " features from " << featureListFile << std::endl;
        }
      }
    }
    if (Globals::VERB > 2) {
      std::cerr << "Partial run: Loading dinosaur features completed" << std::endl;
    }
  }

  /*****************************************************************************
    Step 2: Run MaRaCluster to form MS2 clusters
   ****************************************************************************/

  /* spectrumToPrecursorMap:
       filled by runMaRaCluster();
       consumed by MaRaClusterIO::parseClustersForRTimePairs() */
  std::string maraclusterSubFolder = "maracluster";
  fs::path spectrumToPrecursorFilePath(outputFolder_);
  spectrumToPrecursorFilePath /= maraclusterSubFolder;
  spectrumToPrecursorFilePath /= fnPrefix_ + ".spectrum_to_precursor_map.dat";
  std::string spectrumToPrecursorFile = spectrumToPrecursorFilePath.string();
  std::string clusterFilePath;
  if (!isPartialRun || !partial2MaRaCluster_.empty()) {
    rc = runMaRaCluster(maraclusterSubFolder, fileList, allFeatures, clusterFilePath,
	     spectrumToPrecursorFile, tmpFilePrefixAlign, partial2MaRaCluster_);
    if (rc != EXIT_SUCCESS) return rc;
  }
  
  if (!partial2MaRaCluster_.empty() && partial2MaRaCluster_ != "batch" && partial2MaRaCluster_ != "cluster") {
    if (Globals::VERB > 2) {
      std::cerr << "Partial run MaRaCluster finished" << std::endl;
    }
    return EXIT_SUCCESS;
  }
  
	/*****************************************************************************
	  Step 3: Create minimum spanning tree of alignments
	 ****************************************************************************/
  
  std::vector<std::pair<int, FilePair> > featureAlignmentQueue;
  fs::path featureAlignQueueFile(outputFolder_);
  featureAlignQueueFile /= "maracluster";
  featureAlignQueueFile /= "featureAlignmentQueue.txt";
  
  AlignRetention alignRetention;
  fs::path featureAlignFile(outputFolder_);
  featureAlignFile /= "maracluster";
  featureAlignFile /= "alignRetention.txt";
  
  SpectrumToPrecursorMap spectrumToPrecursorMap(fileList.size());
  
  if (!isPartialRun || (!partial2MaRaCluster_.empty() && !fs::exists(featureAlignFile) && !fs::exists(featureAlignQueueFile))) {
    std::ifstream fileStream(clusterFilePath.c_str(), ios::in);
    
    spectrumToPrecursorMap.deserialize(spectrumToPrecursorFile);
    if (Globals::VERB > 1) {
      std::cerr << "Deserialized spectrum to precursor map" << std::endl;
    }
    
    MaRaClusterIO::parseClustersForRTimePairs(fileStream, fileList, spectrumToPrecursorMap, alignRetention.getRTimePairsRef());

    alignRetention.getAlignModelsTree();
    alignRetention.createMinDepthTree(featureAlignmentQueue);

    if (!partial2MaRaCluster_.empty()) {  // only runs maracluster and outputs the queue to a file
      if (Globals::VERB > 2) {
        std::cerr << "Partial run MaRaCluster: Saving alignments and alignment tree" << std::endl;
      }
      
      rc = saveAlignmentQueue(featureAlignmentQueue, 
                              featureAlignQueueFile.string(), fileList);
      if (rc != EXIT_SUCCESS) return rc;
      
      rc = alignRetention.saveState(featureAlignFile.string());
      if (rc != EXIT_SUCCESS) return rc;
    }
  }
  
  if (!partial2MaRaCluster_.empty()) {
    if (Globals::VERB > 2) {
      std::cerr << "Partial run MaRaCluster finished" << std::endl;
    }
    return EXIT_SUCCESS;
  }
  
  fs::path matchedFeaturesFile(outputFolder_);
  matchedFeaturesFile /= "percolator";
  matchedFeaturesFile /= "matches";
  
  std::string addedFeaturesFile = "";
  if (isPartialRun) {
    std::cerr << "Partial run: Loading alignments and alignment tree" << std::endl;

    // Load state of featureAlignmentQueue
    std::ifstream infile(featureAlignQueueFile.string().c_str(), ios::in);
    int round, prevRound = -1;
    std::string tmp;
    int alignFromIdx, alignToIdx;
    std::vector<bool> featuresAdded(fileList.size());
    size_t originalNumFeaturesFromFile = 0;
    while (infile >> round >> tmp >> tmp >> alignFromIdx >> alignToIdx) {
      if (round != prevRound) {
        prevRound = round;
        originalNumFeaturesFromFile = 0;
      }
      
      FilePair filePair(alignFromIdx, alignToIdx);
      
      fs::path file1(fileList.getFilePath(alignFromIdx));
      fs::path file2(fileList.getFilePath(alignToIdx));
      std::string alignFromFile(file1.filename().string());
      std::string alignToFile(file2.filename().string());

      // Parallel_3 injection. Will only run a pair if it exists in pair/ folder and is as assigned round
      // NOTE: In AlignRetention.cpp, line 197: A run is limited by which round it is, needs to be in order
      if (partial3MatchRound_ == 0 || !partial4Consensus_.empty()) {
        featureAlignmentQueue.push_back(std::make_pair(round, filePair));
      } else if (fs::exists("pair/file1/" + alignFromFile) &&
                 fs::exists("pair/file2/" + alignToFile)) {
        featureAlignmentQueue.push_back(std::make_pair(round, filePair));
        alignRetention.loadState(featureAlignFile.string(), filePair);
        addedFeaturesFile = getAddedFeaturesFN(matchedFeaturesFile, alignFromIdx, alignToIdx);
        if (featuresAdded[alignFromIdx]) allFeatures.at(alignFromIdx).sortByPrecMz();
        if (featuresAdded[alignToIdx]) allFeatures.at(alignToIdx).sortByPrecMz();
        
        if (useTempFiles_) {
          std::string fileName = getFeatureFN(fs::path(tmpFilePrefixAlign) / "features", alignFromIdx);
          std::cerr << "Writing updated features to " << fileName << std::endl;
          bool append = false;
          allFeatures.at(alignFromIdx).saveToFile(fileName, append);
        }
      } else if (round < partial3MatchRound_ &&
          (fs::exists("pair/file1/" + alignToFile) ||
           fs::exists("pair/file2/" + alignToFile))) {
        std::string savedAddedFeaturesFile = getAddedFeaturesFN(matchedFeaturesFile, alignFromIdx, alignToIdx);
        bool withIdxMap = true;
        bool ignoreDuplicates = true;
        bool removePlaceholders = false;
        size_t addedFts = allFeatures.at(alignToIdx).loadFromFile(
            savedAddedFeaturesFile, withIdxMap, ignoreDuplicates, removePlaceholders);
        featuresAdded[alignToIdx] = (addedFts > 0);
        
        if (Globals::VERB > 2) {
          std::cerr << "Read in " << addedFts << " features from " << savedAddedFeaturesFile 
                    << " (total: " << allFeatures.at(alignToIdx).size() << " features)" << std::endl;
        }
        
        // update --use-tmp-files match files for matchFrom run 
        if (fs::exists("pair/file1/" + alignToFile) && useTempFiles_) {
          if (originalNumFeaturesFromFile == 0) {
            originalNumFeaturesFromFile = allFeatures.at(alignToIdx).size() - addedFts;
          }
          
          DinosaurFeatureList addedFeatures;
          
          addedFeatures.loadFromFile(savedAddedFeaturesFile, withIdxMap, 
                                     ignoreDuplicates, removePlaceholders);
          
          updateMatchesTmpFile(fs::path(tmpFilePrefixAlign) / "matches", 
              filePair, allFeatures.at(alignToIdx), addedFeatures,
              originalNumFeaturesFromFile);
        }
      }
    }
    infile.close();
    
    // Load state of alignRetention
    if (partial3MatchRound_ == 0 || (!useTempFiles_  && (partial4Consensus_ == "batch" || partial4Consensus_ == "cluster"))) {
      FilePair alignedFilePair(-1,-1);
      alignRetention.loadState(featureAlignFile.string(), alignedFilePair);
    }
    std::cerr << "Partial run: Loading alignments and alignment tree completed" << std::endl;
  }

	/*****************************************************************************
    Step 4: Match features between runs
   ****************************************************************************/

  fs::path percolatorFolder(rootPath);
  percolatorFolder /= "percolator";
  rc = createDirectory(percolatorFolder);
  if (rc != EXIT_SUCCESS) return rc;

	FeatureAlignment featureAlignment(percolatorFolder.string(), percolatorArgs_,
		      alignPpmTol_, alignRTimeStdevTol_, decoyOffset_, linkPEPThreshold_,
		      linkPEPMbrSearchThreshold_, maxFeatureCandidates_);
  
  /**
    If useTempFiles is set and partial runs are activated, we load the matches 
    from temporary files. Otherwise, we need to get the matches only in "batch"
    and "cluster" modes.   
   **/
  if (!isPartialRun || partial3MatchRound_ >= 0 || (!useTempFiles_ 
      && (partial4Consensus_ == "batch" || partial4Consensus_ == "cluster"))) {
    featureAlignment.matchFeatures(featureAlignmentQueue, fileList, 
        alignRetention, allFeatures, addedFeaturesFile, tmpFilePrefixAlign);
  }

  if (partial3MatchRound_ >= 0) {  // parallel 3 runs dinosaur.
    std::cerr << "Partial run for matching features completed" << std::endl;
    return EXIT_SUCCESS;
  }
  
  /* sort features by index before feature groups processing */
  std::vector<DinosaurFeatureList>::iterator ftListIt;
  for (ftListIt = allFeatures.begin(); ftListIt != allFeatures.end(); ++ftListIt) {
    ftListIt->clearFeatureToIdxMap();
  }

	/*****************************************************************************
	  Step 5: Run MaRaCluster again, now with newly discovered features from
			      targeted Dinosaur runs
	 ****************************************************************************/

	std::string maraclusterSubFolderExtraFeatures = "maracluster_extra_features";
	std::string clusterFilePathExtraFeatures;
	fs::path spectrumToPrecursorFilePathExtraFeatures(outputFolder_);
  spectrumToPrecursorFilePathExtraFeatures /= maraclusterSubFolderExtraFeatures;
  spectrumToPrecursorFilePathExtraFeatures /= fnPrefix_ + ".spectrum_to_precursor_map.dat";
  std::string spectrumToPrecursorFileExtraFeatures = spectrumToPrecursorFilePathExtraFeatures.string();
	rc = runMaRaCluster(maraclusterSubFolderExtraFeatures, fileList, allFeatures,
			   clusterFilePathExtraFeatures, spectrumToPrecursorFileExtraFeatures,
			   tmpFilePrefixAlign, partial4Consensus_);
  if (rc != EXIT_SUCCESS) return rc;
  
  if (!partial4Consensus_.empty() && partial4Consensus_ != "batch" && partial4Consensus_ != "cluster") {
    if (Globals::VERB > 2) {
      std::cerr << "Partial run MaRaCluster finished" << std::endl;
    }
    return EXIT_SUCCESS;
  }
  
  /*****************************************************************************
    Step 6: Apply single linkage clustering to form MS1 feature groups
   ****************************************************************************/

  FeatureGroups featureGroups(maxMissingValues_, intensityScoreThreshold_);
  std::string tmpFilePrefixGroup = "";
  if (useTempFiles_) {
  	fs::path tmpFileFolder(rootPath);
	  tmpFileFolder /= "tmp";
	  tmpFileFolder /= "featureToGroupMaps";
	  rc = createDirectory(tmpFileFolder);
    if (rc != EXIT_SUCCESS) return rc;
    
	  tmpFilePrefixGroup = tmpFileFolder.string() + "/featureToGroupMap.";
  }

	featureGroups.singleLinkClustering(featureAlignmentQueue,
      featureAlignment.getFeatureMatches(), tmpFilePrefixGroup, tmpFilePrefixAlign);
  
	/*****************************************************************************
			Step 7: Keep top 3 consensus spectra per feature group based on
			        intensity score
		****************************************************************************/
  
  if (isPartialRun) {
    spectrumToPrecursorMap.deserialize(spectrumToPrecursorFile);
    if (Globals::VERB > 1) {
      std::cerr << "Deserialized spectrum to precursor map" << std::endl;
    }
  }
  
  std::map<FeatureId, std::vector<int> > featureToSpectrumCluster;
  std::ifstream fileStreamExtraFeatures(clusterFilePathExtraFeatures.c_str(), ios::in);
  MaRaClusterIO::parseClustersForConsensusMerge(fileStreamExtraFeatures,
      fileList, spectrumToPrecursorMap, featureToSpectrumCluster);
  /* this would also link spectra to features found by matches-between-runs,
     but it seems to (slightly) lower specificity
  MaRaClusterIO::parseClustersForConsensusMerge(fileStreamExtraFeatures,
      fileList, spectrumToPrecursorMapExtraFeatures,
      featureToSpectrumCluster);
  */

  if (useTempFiles_) {
    loadAllFeatures(tmpFilePrefixAlign, allFeatures);
  }

  for (ftListIt = allFeatures.begin(); ftListIt != allFeatures.end(); ++ftListIt) {
    ftListIt->buildFeatureIdxToVectorIdxMap();
  }

  featureGroups.filterConsensusFeatures(allFeatures, featureToSpectrumCluster);

  /*****************************************************************************
    Step 8: Write MS1 feature groups to output file
   ****************************************************************************/

  std::vector<DinosaurFeatureList> noFeatures;
  SpectrumToPrecursorMap noMap(fileList.size());
  MaRaClusterAdapter maraclusterAdapter(noFeatures, noMap, outputSpectrumFile_);

  fs::path featureGroupsOutFile(rootPath);
  featureGroupsOutFile /= fnPrefix_ + ".feature_groups.tsv";
  featureGroups.printFeatureGroups(featureGroupsOutFile.string(),
      allFeatures, featureToSpectrumCluster,
      maraclusterAdapter.getSpectrumClusterToConsensusFeatures());
  featureGroups.clear(); // unload feature groups from memory

  if (useTempFiles_) {
    unloadAllFeatures(allFeatures);
  }

  /*****************************************************************************
    Step 9: Create and write consensus spectra
   ****************************************************************************/

  std::string maraclusterSubFolderConsensus = "consensus_spectra";
  std::vector<std::string> maraclusterArgs = maraclusterArgs_;
  fs::path maraclusterFolder(rootPath);
  maraclusterFolder /= maraclusterSubFolderConsensus;
  maraclusterArgs[1] = "consensus";
  maraclusterArgs.push_back("--output-folder");
  maraclusterArgs.push_back(maraclusterFolder.string());
  maraclusterArgs.push_back("--clusterFile");
  maraclusterArgs.push_back(clusterFilePathExtraFeatures);

  maraclusterAdapter.parseOptions(maraclusterArgs);
  maraclusterAdapter.run();

  time_t endTime;
  clock_t endClock = clock();
  time(&endTime);
  double diff_time = difftime(endTime, startTime);

  std::cerr << "Running Quandenser took: "
    << ((double)(endClock - startClock)) / (double)CLOCKS_PER_SEC
    << " cpu seconds or " << diff_time << " seconds wall time" << std::endl;

  return EXIT_SUCCESS;
}

} /* namespace quandenser */
