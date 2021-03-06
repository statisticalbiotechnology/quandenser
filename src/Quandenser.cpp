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

namespace quandenser {

Quandenser::Quandenser() : call_(""), fnPrefix_("Quandenser"), seed_(1u),
    numThreads_(4u), maxMissingValues_(-1), intensityScoreThreshold_(0.5),
    spectrumBatchFileFN_(""), outputFolder_("Quandenser_output"),
    outputSpectrumFile_(""), maraclusterPpmTol_(20.0f),
    alignPpmTol_(20.0f), alignRTimeStdevTol_(10.0f), 
    decoyOffset_(5.0 * 1.000508),
    linkPEPThreshold_(0.25), linkPEPMbrSearchThreshold_(0.05),
    maxFeatureCandidates_(2), useTempFiles_(false) {}

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

void Quandenser::runMaRaCluster(const std::string& maRaClusterSubFolder, 
    const maracluster::SpectrumFileList& fileList, 
    std::vector<DinosaurFeatureList>& allFeatures,
    std::string& clusterFilePath,
    SpectrumToPrecursorMap& spectrumToPrecursorMap,
    const std::string& tmpFilePrefixAlign) {
  std::vector<std::string> maraclusterArgs = maraclusterArgs_;
  boost::filesystem::path maraclusterFolder(outputFolder_);
  maraclusterFolder /= maRaClusterSubFolder;
  maraclusterArgs[1] = "batch";
  maraclusterArgs.push_back("--output-folder");
  maraclusterArgs.push_back(maraclusterFolder.string());
  
  std::string spectrumOutputFile = ""; /* Not needed here */
  MaRaClusterAdapter maraclusterAdapter(allFeatures, spectrumToPrecursorMap, spectrumOutputFile);
  maraclusterAdapter.parseOptions(maraclusterArgs);
  
  boost::filesystem::path spectrumToPrecursorFilePath(maraclusterFolder);
  spectrumToPrecursorFilePath /= fnPrefix_ + ".spectrum_to_precursor_map.dat";
  std::string spectrumToPrecursorFile = spectrumToPrecursorFilePath.string();
  
  if (!boost::filesystem::exists(spectrumToPrecursorFile)) {
    if (useTempFiles_) {
      loadAllFeatures(tmpFilePrefixAlign, allFeatures);
    }
    
    maraclusterAdapter.run();
    
    spectrumToPrecursorMap.serialize(spectrumToPrecursorFile);
    if (Globals::VERB > 1) {
      std::cerr << "Serialized spectrum to precursor map" << std::endl;
    }
    
    if (useTempFiles_) {
      unloadAllFeatures(allFeatures);
    }
  } else {
    spectrumToPrecursorMap.deserialize(spectrumToPrecursorFile);
    if (Globals::VERB > 1) {
      std::cerr << "Deserialized spectrum to precursor map" << std::endl;
    }
  }
  
  clusterFilePath = maraclusterAdapter.getClusterFileName();
}

void Quandenser::loadAllFeatures(const std::string& tmpFilePrefixAlign,
    std::vector<DinosaurFeatureList>& allFeatures) {
  for (size_t i = 0; i < allFeatures.size(); ++i) {
    std::string fileName = tmpFilePrefixAlign + "/features."  + boost::lexical_cast<std::string>(i) + ".dat";
    allFeatures.at(i).loadFromFile(fileName);
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
  
  boost::filesystem::path rootPath(outputFolder_);
  boost::filesystem::create_directories(rootPath, returnedError);
  if (!boost::filesystem::exists(rootPath)) {
    std::cerr << "Error: could not create output directory at " 
        << outputFolder_ << std::endl;
    return EXIT_FAILURE;
  }
  
  /*****************************************************************************  
    Step 0: Read in the input text file with a list of mzML files
   ****************************************************************************/
  
  maracluster::SpectrumFileList fileList;
  fileList.initFromFile(spectrumBatchFileFN_);
  
  if (fileList.size() < 2u) {
    std::cerr << "Error: less than 2 spectrum files were specified, "
        << "Quandenser needs at least two files to perform a "
        << "meaningful alignment." << std::endl;
    return EXIT_FAILURE;
  }
  
  /*****************************************************************************  
    Step 1: Detect MS1 features with Dinosaur
   ****************************************************************************/
  
  boost::filesystem::path dinosaurFolder(rootPath);
  dinosaurFolder /= "dinosaur";
  boost::filesystem::create_directories(dinosaurFolder, returnedError);
  if (!boost::filesystem::exists(dinosaurFolder)) {
    std::cerr << "Error: could not create output directory at " 
        << dinosaurFolder.string() << std::endl;
    return EXIT_FAILURE;
  }
  
  std::vector<boost::filesystem::path> dinosaurFeatureFiles;
  
  std::vector<std::string> files = fileList.getFilePaths();
  for (std::vector<std::string>::const_iterator it = files.begin(); it != files.end(); ++it) {
    if (maracluster::MSFileHandler::getOutputFormat(*it) != "mzml") {
      std::cerr << "Warning: file extension is not .mzML, ignoring file: " << *it << std::endl;
    }
    boost::filesystem::path mzMLFile(*it);
    
    boost::filesystem::path dinosaurFeatureFile(dinosaurFolder.string());
    
    /* TODO: what if we use files with the same filename but in different folders? */
    dinosaurFeatureFile /= mzMLFile.stem().string() + ".features.tsv";
    dinosaurFeatureFiles.push_back(dinosaurFeatureFile);
    if (!boost::filesystem::exists(dinosaurFeatureFile)) {
      if (Globals::VERB > 1) {
        std::cerr << "Processing " << mzMLFile.filename() << " with Dinosaur." << std::endl;
      }
      int rc = DinosaurIO::runDinosaurGlobal(dinosaurFolder.string(), mzMLFile.string());
      if (rc != EXIT_SUCCESS) {
        std::cerr << "Dinosaur failed with exit code " << rc << ". Terminating.." << std::endl;
        return EXIT_FAILURE;
      }
    } else {
      if (Globals::VERB > 1) {
        std::cerr << "Already processed " << mzMLFile.filename() << " with Dinosaur." << std::endl;
      }
    }
  }
  
  std::vector<DinosaurFeatureList> allFeatures;
  size_t fileIdx = 0u;
  for (std::vector<boost::filesystem::path>::const_iterator it = dinosaurFeatureFiles.begin(); it != dinosaurFeatureFiles.end(); ++it) {
    std::ifstream fileStream(it->string().c_str(), ios::in);
    DinosaurFeatureList features;
    DinosaurIO::parseDinosaurFeatureFile(fileStream, fileIdx++, features);
    
    allFeatures.push_back(features);
    if (Globals::VERB > 2) {
      std::cerr << "Read in " << features.size() << " features from " << it->filename() << std::endl;
    }
  }  
  
  std::string tmpFilePrefixAlign = "";
  if (useTempFiles_) {
    boost::filesystem::path tmpFileFolder(rootPath);
    tmpFileFolder /= "tmp";
    tmpFileFolder /= "matchFeatures";
    boost::filesystem::create_directories(tmpFileFolder, returnedError);
    if (!boost::filesystem::exists(tmpFileFolder)) {
      std::cerr << "Error: could not create output directory at " << tmpFileFolder << std::endl;
      return EXIT_FAILURE;
    }
    tmpFilePrefixAlign = tmpFileFolder.string();
    
    for (size_t i = 0; i < allFeatures.size(); ++i) {
      std::string fileName = tmpFilePrefixAlign + "/features."  + boost::lexical_cast<std::string>(i) + ".dat";
      remove(fileName.c_str());
      bool append = false;
      allFeatures.at(i).saveToFile(fileName, append);
      allFeatures.at(i).clear();
    }
  }
  
  /*****************************************************************************  
    Step 2: Run MaRaCluster to form MS2 clusters
   ****************************************************************************/
   
  /* spectrumToPrecursorMap: 
       filled by runMaRaCluster(); 
       consumed by MaRaClusterIO::parseClustersForRTimePairs() */
  SpectrumToPrecursorMap spectrumToPrecursorMap(fileList.size());
  std::string maraclusterSubFolder = "maracluster";
  std::string clusterFilePath;
  runMaRaCluster(maraclusterSubFolder, fileList, allFeatures, clusterFilePath, 
      spectrumToPrecursorMap, tmpFilePrefixAlign);
  
  /*****************************************************************************  
    Step 3: Create minimum spanning tree of alignments
   ****************************************************************************/
   
  AlignRetention alignRetention;
  std::ifstream fileStream(clusterFilePath.c_str(), ios::in);
  MaRaClusterIO::parseClustersForRTimePairs(fileStream, fileList, 
      spectrumToPrecursorMap, alignRetention.getRTimePairsRef());
  
  alignRetention.getAlignModelsTree();
  
  std::vector<std::pair<int, FilePair> > featureAlignmentQueue;
  alignRetention.createMinDepthTree(featureAlignmentQueue);
  
  /*****************************************************************************  
    Step 4: Match features between runs
   ****************************************************************************/
   
  boost::filesystem::path percolatorFolder(rootPath);
  percolatorFolder /= "percolator";
  boost::filesystem::create_directories(percolatorFolder, returnedError);
  if (!boost::filesystem::exists(percolatorFolder)) {
    std::cerr << "Error: could not create output directory at " 
              << percolatorFolder.string() << std::endl;
    return EXIT_FAILURE;
  }
  
  FeatureAlignment featureAlignment(percolatorFolder.string(), percolatorArgs_, 
      alignPpmTol_, alignRTimeStdevTol_, decoyOffset_, linkPEPThreshold_, 
      linkPEPMbrSearchThreshold_, maxFeatureCandidates_);
  featureAlignment.matchFeatures(featureAlignmentQueue, fileList, 
      alignRetention, allFeatures, tmpFilePrefixAlign);
  
  std::vector<DinosaurFeatureList>::iterator ftListIt;
  for (ftListIt = allFeatures.begin(); ftListIt != allFeatures.end(); ++ftListIt) {
    ftListIt->clearFeatureToIdxMap();
  }
  
  /*****************************************************************************  
    Step 5: Apply single linkage clustering to form MS1 feature groups
   ****************************************************************************/
  
  if (maxMissingValues_ < 0) {
    maxMissingValues_ = fileList.size() / 4;
  }
  
  FeatureGroups featureGroups(maxMissingValues_, intensityScoreThreshold_);
  std::string tmpFilePrefixGroup = "";
  if (useTempFiles_) {
    boost::filesystem::path tmpFileFolder(rootPath);
    tmpFileFolder /= "tmp";
    tmpFileFolder /= "featureToGroupMaps";
    boost::filesystem::create_directories(tmpFileFolder, returnedError);
    if (!boost::filesystem::exists(tmpFileFolder)) {
      std::cerr << "Error: could not create output directory at " 
                << tmpFileFolder << std::endl;
      return EXIT_FAILURE;
    }
    tmpFilePrefixGroup = tmpFileFolder.string() + "/featureToGroupMap.";
  }
  featureGroups.singleLinkClustering(featureAlignmentQueue, 
      featureAlignment.getFeatureMatches(), tmpFilePrefixGroup, tmpFilePrefixAlign);
  
  /*****************************************************************************  
    Step 6: Run MaRaCluster again, now with newly discovered features from
            targeted Dinosaur runs
   ****************************************************************************/
  
  std::string maraclusterSubFolderExtraFeatures = "maracluster_extra_features";
  std::string clusterFilePathExtraFeatures;
  SpectrumToPrecursorMap spectrumToPrecursorMapExtraFeatures(fileList.size());
  runMaRaCluster(maraclusterSubFolderExtraFeatures, fileList, allFeatures, 
      clusterFilePathExtraFeatures, spectrumToPrecursorMapExtraFeatures, 
      tmpFilePrefixAlign);  
      
  /*****************************************************************************  
    Step 7: Keep top 3 consensus spectra per feature group based on 
            intensity score 
   ****************************************************************************/
  
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
    ftListIt->sortByFeatureIdx();
  }
  
  featureGroups.filterConsensusFeatures(allFeatures, featureToSpectrumCluster);
  
  /*****************************************************************************  
    Step 8: Write MS1 feature groups to output file
   ****************************************************************************/
  
  std::vector<DinosaurFeatureList> noFeatures;
  SpectrumToPrecursorMap noMap(fileList.size());
  MaRaClusterAdapter maraclusterAdapter(noFeatures, noMap, outputSpectrumFile_);
  
  boost::filesystem::path featureGroupsOutFile(rootPath);
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
  boost::filesystem::path maraclusterFolder(rootPath);
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
