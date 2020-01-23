/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

#include "Option.h"

namespace quandenser {

const std::string Option::NO_SHORT_OPT = "NO_SHORT_OPT_CONSTANT";
const std::string Option::EXPERIMENTAL_FEATURE = "EXPERIMENTAL_FEATURE_CONSTANT";

template<class T>
bool from_string(T& t, const std::string& s) {
  std::istringstream iss(s);
  return !(iss >> t).fail();
}

void searchandreplace(std::string& source, const std::string& find,
                      const std::string& replace) {
  size_t j;
  for (; (j = source.find(find)) != std::string::npos;) {
    source.replace(j, find.length(), replace);
  }
}

Option::Option(std::string shrt, std::string lng, std::string nm, std::string hlp,
               std::string hlpType, OptionOption typ, std::string dfl) {
  type = typ;
  shortOpt = shrt;
  longOpt = lng;
  help = hlp;
  helpType = hlpType;
  name = nm;
  deflt = dfl;
}

Option::~Option() {
}

bool Option::operator ==(const std::string& option) {
  return ((shortOpt != Option::NO_SHORT_OPT && 
           shortOpt != Option::EXPERIMENTAL_FEATURE && 
           shortOpt == option) 
          || longOpt == option);
}

CommandLineParser::CommandLineParser(std::string usage, std::string tail) {
  header = usage;
  endnote = tail;
  optMaxLen = 0;
  defineOption("h", "help", "Display this message");
}

CommandLineParser::~CommandLineParser() {
}

double CommandLineParser::getDouble(std::string dest, double lower,
                                    double upper) {
  double val;
  from_string<double> (val, options[dest]);
  if (!from_string<double> (val, options[dest]) || (val < lower || val > upper)) 
  {
    std::ostringstream temp;
    temp << "-" << dest << " option requires a float between " << lower
        << " and " << upper << std::endl;
    throw MyException(temp.str());
  }
  return val;
}

int CommandLineParser::getInt(std::string dest, int lower, int upper) {
  int val;
  if (!from_string<int> (val, options[dest]) || val < lower || val > upper) 
  {
    std::ostringstream temp;
    temp << "-" << dest << " option requires an integer between " << lower
        << " and " << upper << std::endl;
    throw MyException(temp.str());
  }
  return val;
}

void CommandLineParser::defineOption(std::string shortOpt, std::string longOpt,
                                     std::string help, std::string helpType,
                                     OptionOption typ, std::string dfault) {
  //NOTE brute force to check if the option is already defined
  for(std::vector<Option>::const_iterator it = opts.begin();
      it != opts.end(); it++) {
	  if((shortOpt != Option::NO_SHORT_OPT && 
	      shortOpt != Option::EXPERIMENTAL_FEATURE && 
	      (*it).shortOpt == shortOpt) 
	     || (*it).longOpt == longOpt) {
	    std::ostringstream temp;
	    temp << "ERROR : option " << shortOpt << "," << longOpt << " is already defined " << std::endl;
	    throw MyException(temp.str());
	  }
  }
  
  opts.insert(opts.begin(), Option((shortOpt == Option::NO_SHORT_OPT || shortOpt == Option::EXPERIMENTAL_FEATURE) ? shortOpt : "-" + shortOpt,
                                   "--" + longOpt,
                                   longOpt,
                                   help,
                                   helpType,
                                   typ,
                                   dfault));
  if (longOpt.length() + helpType.length() > optMaxLen) {
    optMaxLen = longOpt.length() + helpType.length();
  }
}

void CommandLineParser::parseArgs(int argc, char** argv) {
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      findOption(argv, i);
    } else {
      arguments.insert(arguments.end(), argv[i]);
    }
  }
}

void CommandLineParser::error(std::string msg) {
  std::ostringstream temp;
  temp << header << std::endl << msg << std::endl;
  throw MyException(temp.str());
}

void CommandLineParser::help() {
  std::string::size_type descLen = optMaxLen + 8;
  std::string::size_type helpLen = lineLen - descLen;
  std::cerr << header << std::endl << "Options:" << std::endl;
  for (size_t i = opts.size(); i--;) {
    std::string::size_type j = 0;
    if (opts[i].shortOpt != Option::NO_SHORT_OPT && opts[i].shortOpt != Option::EXPERIMENTAL_FEATURE) {
      std::cerr << " " << opts[i].shortOpt;
      if (opts[i].helpType.length() > 0) {
        std::cerr << " <" << opts[i].helpType << ">";
      }
    } else if (opts[i].shortOpt == Option::EXPERIMENTAL_FEATURE) {
      std::cerr << "[EXPERIMENTAL FEATURE]";
    }
    std::cerr << std::endl;
    std::string desc = " " + opts[i].longOpt;
    if (opts[i].helpType.length() > 0) {
      desc += " <" + opts[i].helpType + ">";
    }
    while (j < opts[i].help.length()) {
      std::cerr.width(descLen);
      std::cerr << std::left << desc;
      desc = " ";
      std::cerr.width(0);
      std::string::size_type l = helpLen;
      if (j + l < opts[i].help.length()) {
        std::string::size_type p = opts[i].help.rfind(' ', j + l);
        if (p != std::string::npos && p > j) {
          l = p - j + 1;
        }
      }
      std::cerr << opts[i].help.substr(j, l) << std::endl;
      j += l;
    }
  }
  std::cerr << std::endl << endnote << std::endl;
  exit(0);
}

void CommandLineParser::htmlHelp() {
  std::cerr << "<html><title>Title</title><body><blockquote>" << std::endl;
  std::cerr
      << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\" />"
      << std::endl;
  std::string htmlHeader = header;
  searchandreplace(htmlHeader, "\n", "<br/>");
  std::cerr << htmlHeader << std::endl << "Options:" << std::endl;
  std::cerr << "<table border=0>" << std::endl;
  for (size_t i = opts.size(); i--;) {
    std::cerr << "<tr><td><code>";
    if (opts[i].shortOpt != Option::NO_SHORT_OPT && opts[i].shortOpt != Option::EXPERIMENTAL_FEATURE) {
      std::cerr << opts[i].shortOpt;
      if (opts[i].helpType.length() > 0) {
        std::cerr << " &lt;" << opts[i].helpType << "&gt;";
      }
      std::cerr << "</code>, ";
    } else if (opts[i].shortOpt == Option::EXPERIMENTAL_FEATURE) {
      std::cerr << "[EXPERIMENTAL FEATURE]";
    }
    std::cerr << "<code>";
    std::cerr << " " + opts[i].longOpt;
    if (opts[i].helpType.length() > 0) {
      std::cerr << " &lt;" << opts[i].helpType << "&gt;";
    }
    std::cerr << "</code></td>" << std::endl;
    std::cerr << "<td>" << opts[i].help << "</td></tr>" << std::endl;
  }
  std::cerr << "</table>" << std::endl;
  std::string htmlEnd = endnote;
  searchandreplace(htmlEnd, "\n", "<br>");
  std::cerr << "<br/>" << std::endl << htmlEnd << "<br/>" << std::endl;
  std::cerr << "</blockquote></body></html>" << std::endl;
  exit(0);
}

void CommandLineParser::findOption(char** argv, int& index) {
  if ((std::string)argv[index] == "-html" || (std::string)argv[index] == "--html") {
    htmlHelp();
  }
  if ((std::string)argv[index] == "-h" || (std::string)argv[index] == "--help") {
    help();
  }
  std::string optstr = (std::string)argv[index];
  std::string valstr("");
  std::string::size_type eqsign = optstr.find('=');
  if (eqsign != std::string::npos) {
    valstr = optstr.substr(eqsign + 1);
    optstr = optstr.substr(0, eqsign);
  }
  for (size_t i = 0; i < opts.size(); i++) {
    if (opts[i] == optstr) {
      switch (opts[i].type) {
        case FALSE_IF_SET:
          options[opts[i].name] = "0";
          break;
        case TRUE_IF_SET:
          options[opts[i].name] = "1";
          break;
        case VALUE:
          if (valstr.length() > 0) {
            options[opts[i].name] = valstr;
          } else {
            options[opts[i].name] = argv[index + 1];
            index++;
          }
          break;
        case MAYBE:
          if (valstr.length() > 0) {
            options[opts[i].name] = valstr;
          } else if (argv[index + 1][0] != '-') {
            options[opts[i].name] = argv[index + 1];
            index++;
          } else {
            options[opts[i].name] = opts[i].deflt;
          }
          break;
        default:
          break;
      };
      return;
    }
  }
  error("ERROR: the option " + optstr + " is invalid.\n" +
        "Please run \"command --help.\"");
}

} /* namespace quandenser */
