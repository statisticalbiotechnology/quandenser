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

#ifndef QUANDENSER_SPLINEREGRESSION_H_
#define QUANDENSER_SPLINEREGRESSION_H_

#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <sstream>

#include "percolator/src/BaseSpline.h"

#include "Globals.h"

namespace quandenser {

class SplineRegression : public BaseSpline {
  public:
    SplineRegression(){};
    virtual ~SplineRegression(){};
    void predict(const std::vector<double>& x, std::vector<double>& predict) {
      return BaseSpline::predict(x, predict);
    }
    void setData(const std::vector<double>& xx, const std::vector<double>& yy) {
      x = xx;
      y = yy;
    }
    std::pair<std::vector<double>, std::vector<double> > getData() {
      std::pair<std::vector<double>, std::vector<double> > vector_pair = std::make_pair(x, y);
      return vector_pair;
    }
    inline void setRmse(const float rmse) { rmse_ = rmse; }
    inline float getRmse() const { return rmse_; }
  protected:
    virtual void calcPZW();
    virtual void initg();
    virtual void limitg();
    virtual void limitgamma();
    std::vector<double> y;
    static const double gRange;

    float rmse_;
};

} /* namespace quandenser */

#endif /* QUANDENSER_SPLINEREGRESSION_H_ */
