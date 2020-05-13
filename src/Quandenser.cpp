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

namespace quandenser {

Quandenser::Quandenser() : call_(""), fnPrefix_("Quandenser"), seed_(1u),
    maxMissingValues_(-1), intensityScoreThreshold_(0.5),
    spectrumBatchFileFN_(""), outputFolder_("Quandenser_output"),
    outputSpectrumFile_(""), maraclusterPpmTol_(20.0f),
    alignPpmTol_(20.0f), alignRTimeStdevTol_(10.0f),
    decoyOffset_(5.0 * 1.000508),
    linkPEPThreshold_(0.25), linkPEPMbrSearchThreshold_(0.05),
    maxFeatureCandidates_(2),  useTempFiles_(false),
    parallel_1_(0), parallel_2_(0), parallel_3_(0), parallel_4_(0) {}

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
      "Maximum number of threads available for Dinosaur (default: 4).",
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
      "parallel-1",
      "Parallel quandenser, stop1",
      "int");
  cmd.defineOption("X",
      "parallel-2",
      "Parallel quandenser, stop2",
      "int");
  cmd.defineOption("Y",
      "parallel-3",
      "Parallel quandenser, stop3",
      "int");
  cmd.defineOption("Z",
      "parallel-4",
      "Parallel quandenser, stop4",
      "int");
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

  if (cmd.optionSet("dinosaur-memory")) {
    DinosaurIO::setJavaMemory(cmd.options["dinosaur-memory"]);
  }

  if (cmd.optionSet("dinosaur-threads")) {
    DinosaurIO::setJavaNumThreads(cmd.getInt("dinosaur-threads", 1, 100));
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

  if (cmd.optionSet("parallel-1")) {
    parallel_1_ = 1;
  }

  if (cmd.optionSet("parallel-2")) {
    parallel_2_ = 1;
  }

  if (cmd.optionSet("parallel-3")) {
    parallel_3_ = cmd.getInt("parallel-3", 0, 99999);  // Maximum depth is 99999
  }

  if (cmd.optionSet("parallel-4")) {
    parallel_4_ = 1;
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

    if (useTempFiles_) {
      unloadAllFeatures(allFeatures);
    }

    spectrumToPrecursorMap.serialize(spectrumToPrecursorFile);
    if (Globals::VERB > 1) {
      std::cerr << "Serialized spectrum to precursor map" << std::endl;
    }
  } else {
    spectrumToPrecursorMap.deserialize(spectrumToPrecursorFile);
    if (Globals::VERB > 1) {
      std::cerr << "Deserialized spectrum to precursor map" << std::endl;
    }
  }

  clusterFilePath = maraclusterAdapter.getClusterFileName();
}

std::string Quandenser::getFeatureFN(
    const boost::filesystem::path& featureOutFile, size_t fileIdx) {
  return featureOutFile.string() + "." +
    boost::lexical_cast<std::string>(fileIdx) + ".dat";
}

std::string Quandenser::getFeatureFN(
    const boost::filesystem::path& featureOutFile,
    size_t fileIdx1, size_t fileIdx2) {
  return featureOutFile.string() + "." +
    boost::lexical_cast<std::string>(fileIdx1) + "_" +
    boost::lexical_cast<std::string>(fileIdx2) + ".dat";
}

void Quandenser::loadAllFeatures(const std::string& tmpFilePrefixAlign,
    std::vector<DinosaurFeatureList>& allFeatures) {
  for (size_t i = 0; i < allFeatures.size(); ++i) {
    std::string fileName = tmpFilePrefixAlign + "/features."  + boost::lexical_cast<std::string>(i) + ".dat";
    allFeatures.at(i).loadFromFile(fileName);
  }
}

void Quandenser::unloadAllFeatures(std::vector<DinosaurFeatureList>& allFeatures) {
  for (size_t i = 0; i < allFeatures.size(); ++i) {
    allFeatures.at(i).clear();
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

  if (fileList.size() < 2u && !parallel_1_) {
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
  if (parallel_1_) {
    std::cout << "Parallel stop 1 reached" << std::endl;
    return EXIT_SUCCESS;
  }

  boost::filesystem::path featureOutFile(outputFolder_);
  featureOutFile /= "dinosaur/features";

  std::vector<DinosaurFeatureList> allFeatures;
  if (!parallel_3_ && !parallel_4_) {
    size_t fileIdx = 0u;
    for (std::vector<boost::filesystem::path>::const_iterator it = dinosaurFeatureFiles.begin(); it != dinosaurFeatureFiles.end(); ++it) {
      std::ifstream fileStream(it->string().c_str(), ios::in);
      DinosaurFeatureList features;
      DinosaurIO::parseDinosaurFeatureFile(fileStream, fileIdx++, features);

      features.sortByPrecMz();

      allFeatures.push_back(features);
      if (Globals::VERB > 2) {
        std::cerr << "Read in " << features.size() << " features from " << it->filename() << std::endl;
      }
    }

    if (parallel_2_) {
      std::cout << "Parallel 2: Saving dinosaur features" << std::endl;
      for (size_t allFeaturesIdx = 0u; allFeaturesIdx < allFeatures.size(); allFeaturesIdx++) {
        DinosaurFeatureList ftList = allFeatures.at(allFeaturesIdx);
        std::string featureListFile = getFeatureFN(featureOutFile, allFeaturesIdx);
        bool append = false;
        ftList.saveToFile(featureListFile, append);
      }
      std::cout << "Parallel 2: Save completed" << std::endl;
    }
  } else {
    std::cout << "Parallel 3/4: Loading dinosaur features from Parallel-2" << std::endl;
    size_t fileIdx = 0u;
    for (std::vector<boost::filesystem::path>::const_iterator it = dinosaurFeatureFiles.begin(); it != dinosaurFeatureFiles.end(); ++it) {
      DinosaurFeatureList ftList;
      allFeatures.push_back(ftList);
      boost::filesystem::path ftFilePath(fileList.getFilePath(fileIdx));
      std::string fn(ftFilePath.filename().string());
      // Check if the pair assigned exists in the pair folder
      if (boost::filesystem::exists("pair/file1/" + fn) ||
          boost::filesystem::exists("pair/file2/" + fn) ||
          parallel_4_) {
        std::string featureListFile = getFeatureFN(featureOutFile, fileIdx);
        size_t ftsAdded = allFeatures.at(fileIdx).loadFromFile(featureListFile);
        std::cerr << "Read in " << allFeatures.at(fileIdx).size() << " features from " << featureListFile << std::endl;
      }
      ++fileIdx;
    }
    std::cout << "Parallel 3/4: Loading completed" << std::endl;
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
  boost::filesystem::path featureAlignFile(outputFolder_);
  featureAlignFile /= "maracluster/featureAlignmentQueue.txt";
  std::vector<std::pair<int, FilePair> > featureAlignmentQueue;
  std::string addedFeaturesFile = "";
  if (!parallel_3_ && !parallel_4_) {
    std::ifstream fileStream(clusterFilePath.c_str(), ios::in);
    MaRaClusterIO::parseClustersForRTimePairs(fileStream, fileList, spectrumToPrecursorMap, alignRetention.getRTimePairsRef());

    alignRetention.getAlignModelsTree();
    alignRetention.createMinDepthTree(featureAlignmentQueue);

    if (parallel_2_) {  // parallel 2 only runs maracluster and outputs the queue to a file
      std::cout << "Parallel 2: Saving alignment" << std::endl;
      // Save featureAlignmentQueue vector, used by Nextflow and quandenser
      std::ofstream outfile(featureAlignFile.string().c_str(), ios::out);
      std::vector<std::pair<int, FilePair> >::const_iterator filePairIt;
      for (filePairIt = featureAlignmentQueue.begin();
            filePairIt != featureAlignmentQueue.end(); ++filePairIt) {
        FilePair filePair = filePairIt->second;
        int round = filePairIt->first;
        boost::filesystem::path file1(fileList.getFilePath(filePair.fileIdx1));
        boost::filesystem::path file2(fileList.getFilePath(filePair.fileIdx2));
        // Serialized output stream
        outfile << round << '\t' << file1.string() << '\t' << file2.string()
                << '\t' << filePair.fileIdx1 << '\t' << filePair.fileIdx2 << std::endl;
      }
      outfile.close();

      // Save alignRetention state
      alignRetention.SaveState();
      std::cout << "Parallel 2: Saving completed" << std::endl;

      std::cout << "Parallel stop 2 reached" << std::endl;
      return EXIT_SUCCESS;
    }
  } else {
    std::cout << "Parallel 3/4: Loading alignment from Parallel 2" << std::endl;

    // Load state of featureAlignmentQueue
    std::ifstream infile(featureAlignFile.string().c_str(), ios::in);
    int round;
    std::string tmp1, tmp2;
    int fileidx1, fileidx2;
    int loop_counter = 0;  // Used for parallelization
    std::vector<bool> featuresAdded(fileList.size());
    while (infile >> round >> tmp1 >> tmp2 >> fileidx1 >> fileidx2) {
      FilePair filePair(fileidx1, fileidx2);

      boost::filesystem::path file1(fileList.getFilePath(fileidx1));
      boost::filesystem::path file2(fileList.getFilePath(fileidx2));
      std::string fn1(file1.filename().string());
      std::string fn2(file2.filename().string());

      // Parallel_3 injection. Will only run a pair if it exists in pair/ folder and is as assigned round
      // NOTE: In AlignRetention.cpp, line 197: A run is limited by which round it is, needs to be in order
      if (parallel_4_) {
        featureAlignmentQueue.push_back(std::make_pair(round, filePair));
      } else if (boost::filesystem::exists("pair/file1/" + fn1) &&
                 boost::filesystem::exists("pair/file2/" + fn2)) {
        featureAlignmentQueue.push_back(std::make_pair(round, filePair));
        addedFeaturesFile = getFeatureFN(featureOutFile, fileidx1, fileidx2);
        if (featuresAdded[fileidx1]) allFeatures.at(fileidx1).sortByPrecMz();
        if (featuresAdded[fileidx2]) allFeatures.at(fileidx2).sortByPrecMz();
      } else if (loop_counter < parallel_3_ - 1 &&  // -1 because parallel_3_ is added +1 from actual depth
          (boost::filesystem::exists("pair/file1/" + fn2) ||
           boost::filesystem::exists("pair/file2/" + fn2))) {
        std::string savedAddedFeaturesFile = getFeatureFN(featureOutFile, fileidx1, fileidx2);
        size_t addedFts = allFeatures.at(fileidx2).loadFromFile(savedAddedFeaturesFile);
        std::cerr << "Read in " << addedFts << " features from " << savedAddedFeaturesFile << std::endl;
        featuresAdded[fileidx2] = true;
      }
      loop_counter++;
    }
    infile.close();

    // Load state of alignRetention
    alignRetention.LoadState();
    std::cout << "Parallel 3/4: Loading completed" << std::endl;
  }

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

  std::string percolatorOutputFileBaseFN = percolatorFolder.string();

  if (maxMissingValues_ < 0) {
    maxMissingValues_ = fileList.size() / 4;
  }

	FeatureAlignment featureAlignment(percolatorFolder.string(), percolatorArgs_,
		      alignPpmTol_, alignRTimeStdevTol_, decoyOffset_, linkPEPThreshold_,
		      linkPEPMbrSearchThreshold_, maxFeatureCandidates_);

  featureAlignment.matchFeatures(featureAlignmentQueue, fileList, alignRetention, allFeatures, addedFeaturesFile, tmpFilePrefixAlign);

  if (parallel_3_) {  // parallel 3 runs dinosaur.
    std::cout << "Parallel stop 3 reached" << std::endl;
    return EXIT_SUCCESS;
  }

  /* new maracluster run with newly added features */
  std::string maraclusterSubFolderExtraFeatures = "maracluster_extra_features";
  std::string clusterFilePathExtraFeatures;
  SpectrumToPrecursorMap spectrumToPrecursorMapExtraFeatures(fileList.size());
  runMaRaCluster(maraclusterSubFolderExtraFeatures, fileList, allFeatures, clusterFilePathExtraFeatures, spectrumToPrecursorMapExtraFeatures);
  if (parallel_4_) {
    //std::cout << "Parallel stop 4 reached" << std::endl;
    //return EXIT_SUCCESS;
  }

  /* sort features by index before feature groups processing */
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
