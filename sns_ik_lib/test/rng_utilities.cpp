/*
 *    Copyright 2018 Rethink Robotics
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

#include "rng_utilities.hpp"

#include <random>
#include <ros/console.h>

namespace sns_ik {
namespace rng_util {

/*************************************************************************************************/

void setRngSeed(int seedDouble, int seedInt) {
  getRngDouble(seedDouble);
  getRngInt(seedInt);
}

/*************************************************************************************************/

double getRngDouble(int seed, double low, double upp) {
  static std::mt19937 gen;  // mersenne-twister pseudo-random number generator
  if (upp <= low) { return low; } // catch edge case
  if (seed > 0) { gen.seed(seed); }  // set the seed so that we get repeatable outputs for the test
  static std::uniform_real_distribution<> dis(0.0, 1.0);
  return low + (upp - low) * dis(gen);  // sample the RNG and then map to the range [0, 1]
}

/*************************************************************************************************/

int getRngInt(int seed, int low, int upp) {
  if (upp <= low) { return low; } // catch edge case
  if (seed != 0) { srand(seed); }  // set the seed so that we get repeatable outputs
  return rand() % (upp-low) + low;  // not a true uniform distribution, but close enough
}

/*************************************************************************************************/

Eigen::MatrixXd getRngMatrixXd(int seed, int nRows, int nCols, double low, double upp)
{
  Eigen::MatrixXd data(nRows, nCols);
  if (nRows < 1 || nCols < 1) { return data; }
  for (int iRow = 0; iRow < nRows; iRow++) {
    for (int iCol = 0; iCol < nCols; iCol++) {
      if (iRow == 0 && iCol == 0) {  // Set the seed in the RNG on first call
       data(iRow, iCol) = getRngDouble(seed, low, upp);
      } else {  // Do not update the seed
       data(iRow, iCol) = getRngDouble(0, low, upp);  // 0 == do not update seed
      }
    }
  }
  return data;
}

/*************************************************************************************************/

Eigen::MatrixXd getRngMatrixXdRanked(int seed, int nRows, int nCols, int nRank)
{
  nRows = std::max(nRows, 1);
  nCols = std::max(nCols, 1);
  nRank = std::max(nRank, 1);
  nRank = std::min({nRows, nCols, nRank});
  Eigen::MatrixXd A = getRngMatrixXd(seed + 67394, nRows, nRank, -1.0, 1.0);
  Eigen::MatrixXd B = getRngMatrixXd(seed + 78895, nRank, nCols, -1.0, 1.0);
  Eigen::MatrixXd X = A * B;
  return X;
}

/*************************************************************************************************/

Eigen::ArrayXd getRngArrBndXd(int seed, const Eigen::ArrayXd& low, const Eigen::ArrayXd& upp)
{
  int n = std::min(low.size(), upp.size());
  Eigen::ArrayXd data(n);
  for (int i = 0; i < n; i++) {
    data(i) = getRngDouble(seed, low(i), upp(i));
    seed = 0;  // automatically update RNG after first call
  }
  return data;
}

/*************************************************************************************************/

KDL::JntArray getRngBoundedJoints(int seed, const KDL::JntArray& qLow, const KDL::JntArray& qUpp)
{
  if (qUpp.rows() != qLow.rows()) { ROS_ERROR("Invalid input!");  return qLow; }
  int nJnt = qLow.rows();
  KDL::JntArray q(nJnt);
  for (int iJnt = 0; iJnt < nJnt; iJnt++) {
    q(iJnt) = getRngDouble(seed, qLow(iJnt), qUpp(iJnt));
    seed = 0;  // let the RNG update the sequence automatically after the first call
  }
  return q;
}

/*************************************************************************************************/

KDL::JntArray getNearbyJoints(int seed, const KDL::JntArray& qNom, double delta,
                              const KDL::JntArray& qLow, const KDL::JntArray& qUpp)
{
  if (qLow.rows() != qNom.rows()) { ROS_ERROR("Invalid input!");  return qNom; }
  if (qUpp.rows() != qNom.rows()) { ROS_ERROR("Invalid input!");  return qNom; }
  int nJnt = qLow.rows();
  KDL::JntArray q(nJnt);
  for (int iJnt = 0; iJnt < nJnt; iJnt++) {
    double low = std::max(qNom(iJnt) - delta, qLow(iJnt));
    double upp = std::min(qNom(iJnt) + delta, qUpp(iJnt));
    q(iJnt) = getRngDouble(seed, low, upp);
    seed = 0;  // let the RNG update the sequence automatically after the first call
  }
  return q;
}

/*************************************************************************************************/

}  // namespace rng_util
}  // namespace sns_ik
