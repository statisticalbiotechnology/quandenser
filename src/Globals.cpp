/******************************************************************************  
  Copyright 2015 Matthew The <matthew.the@scilifelab.se>
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
 
#include "Globals.h"
#include "MyException.h"
#include "Version.h"

namespace quandenser {

unsigned int Globals::VERB = 3;

bool Globals::fileExists(const std::string& fileName) {
  std::ifstream infile(fileName.c_str());
  return infile.good();
}

bool Globals::fileIsEmpty(const std::string& fileName) {
  std::ifstream in(fileName.c_str(), std::ios::ate | std::ios::binary);
  if (in.is_open()) {
    return static_cast<long long>(in.tellg()) == 0ll;
  } else {
    return true;
  }
}

void Globals::reportProgress(time_t& startTime, clock_t& startClock,
    size_t currentIt, size_t totalIt) {
  time_t elapsedTime;
  time(&elapsedTime);
  clock_t elapsedClock = clock();
  
  double diff = difftime(elapsedTime, startTime);
  
  unsigned int timeElapsedMin = static_cast<unsigned int>(diff/60);
  unsigned int timeElapsedSecMod = 
      static_cast<unsigned int>(diff - timeElapsedMin * 60);
  
  double elapsedCpuTime = (elapsedClock - startClock) / (double)CLOCKS_PER_SEC;
  std::cerr << "  Elapsed time: " << elapsedCpuTime << " cpu seconds " <<
               "or " << timeElapsedMin << " min " << timeElapsedSecMod << 
               " sec wall time." << std::endl;
  
  double timeLeftSec = (diff / (currentIt+1)) * (totalIt - (currentIt+1));
  unsigned int timeLeftMin = static_cast<unsigned int>(timeLeftSec/60);
  unsigned int timeLeftSecMod = 
      static_cast<unsigned int>(timeLeftSec - timeLeftMin * 60);
  std::cerr << "  Estimated time remaining: " << timeLeftMin << " min " <<
               timeLeftSecMod << " sec wall time." << std::endl;
}


#if defined (__WIN32__) || defined (__MINGW__) || defined (MINGW) || defined (_WIN32)
#include <windows.h>
#include <tchar.h>
#endif

const std::string Globals::getJarPath() {
  std::string out = JAR_PATH;
#if defined (__WIN32__) || defined (__MINGW__) || defined (MINGW) || defined (_WIN32)
  std::wstring keyName = L"Software\\Quandenser\\quandenser-";
  keyName += LVERSION_NAME;
  HKEY hKey;
  RegOpenKeyExW(HKEY_LOCAL_MACHINE, keyName.c_str(), 0, KEY_READ, &hKey);

  WCHAR szBuffer[512];
  DWORD dwBufferSize = sizeof(szBuffer);
  ULONG nError;
  std::wstring strValueName = L"";
  nError = RegQueryValueExW(hKey, strValueName.c_str(), 0, NULL, (LPBYTE)szBuffer, &dwBufferSize);
  if (ERROR_SUCCESS == nError)
  {
    char szcBuffer[512];
    char DefChar = ' ';
    WideCharToMultiByte(CP_ACP,0,szBuffer,-1, szcBuffer,512,&DefChar, NULL);
    out = szcBuffer;
    out += "\\";
    out += JAR_PATH;
  } else {
    
    keyName = L"Software\\Wow6432Node\\Quandenser\\quandenser-";
    keyName += LVERSION_NAME;
    RegOpenKeyExW(HKEY_LOCAL_MACHINE, keyName.c_str(), 0, KEY_READ, &hKey);

    WCHAR szBuffer[512];
    DWORD dwBufferSize = sizeof(szBuffer);
    ULONG nError;
    std::wstring strValueName = L"";
    nError = RegQueryValueExW(hKey, strValueName.c_str(), 0, NULL, (LPBYTE)szBuffer, &dwBufferSize);
    if (ERROR_SUCCESS == nError)
    {
      char szcBuffer[512];
      char DefChar = ' ';
      WideCharToMultiByte(CP_ACP,0,szBuffer,-1, szcBuffer,512,&DefChar, NULL);
      out = szcBuffer;
      out += "\\";
      out += JAR_PATH;
    }
  }
#endif
  return out;  
}


} /* namespace quandenser */
