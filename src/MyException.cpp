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

#include "MyException.h"

namespace quandenser {

MyException::MyException(const std::string& ss) : msg(ss) 
{
  
}

MyException::MyException(const std::ostream& ss)
{
   std::ostringstream x;
   x << ss.rdbuf();
   msg = std::string(x.str());
}


MyException::~MyException() throw()
{

}

const char* MyException::what() const throw()
{
  return msg.c_str(); 
}

} /* namespace quandenser */
