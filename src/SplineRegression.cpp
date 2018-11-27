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

#include "SplineRegression.h"

namespace quandenser {

void SplineRegression::limitg() {
}

void SplineRegression::limitgamma() {
}

void SplineRegression::calcPZW() {
  for (int ix = z.numberEntries(); ix--;) {
    assert(isfinite(g[ix]));
    double epsilon = 1e-4;
    w.packedReplace(ix, 1.0 / max(std::abs(y[ix] - g[ix]), epsilon));
    assert(isfinite(w[ix]));
    z.packedReplace(ix, y[ix]);
    assert(isfinite(z[ix]));
  }
}

void SplineRegression::initg() {
  BaseSpline::initg();
  int n = x.size();
  gnew = PackedVector(n);
  for (int ix = g.size(); ix--;) {
    gnew.packedReplace(ix, y[g.size() / 2]);
    assert(isfinite(gnew[ix]));
  }
}


} /* namespace quandenser */
