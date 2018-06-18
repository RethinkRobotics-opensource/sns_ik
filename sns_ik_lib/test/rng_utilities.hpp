/** @file rng_utilities.hpp
 *
 * @brief a collection of functions for generating repeatable pseudo-random test data
 *
 * @author Matthew Kelly
 *
 * This file provides a set of functions that are used to generate repeatable pseudo-random
 * data for unit tests. All functions provide the caller with direct control over the seed that
 * is used by the random number generator. Passing a seed of zero is used to indicate that the
 * seed should not be reset, instead letting the random generator progress to the next value in
 * the sequence. This is useful for generating numbers in a loop, such as for a matrix.
 *
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

#ifndef SNS_IK_RNG_UTILITIES_H_
#define SNS_IK_RNG_UTILITIES_H_

#include <Eigen/Dense>
#include <kdl/jntarray.hpp>

namespace sns_ik {
namespace rng_util {

/*************************************************************************************************/

/*
 * Set the seeds that are used for both random number generators in this library.
 * @param seedDouble: the seed that is used for generating double precision floating point data
 * @param seedInt: the seed that is used for generating int data
 */
void setRngSeed(int seedDouble, int seedInt);

/*************************************************************************************************/

/*
 * Generate a pseudo-random value on the range [low, upp], given a starting seed.
 * Given the same input seed, it will return the same output.
 * The distribution is uniform between the specified bounds.
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param[opt] low: lower bound on the output range  (default 0.0)
 * @param[opt] upp: upper bound on the output range  (default 1.0)
 * @return: a pseudorandom value between low and upp  (return low if upp <= low)
 */
double getRngDouble(int seed = 0, double low = 0.0, double upp = 1.0);

/*************************************************************************************************/

/*
 * Generate a pseudo-random value on the range [low, upp], given a starting seed.
 * Given the same input seed, it will return the same output
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param[opt] low: lower bound on the output range  (default 0)
 * @param[opt] upp: upper bound on the output range  (default 9)
 * @return: a pseudorandom value between low and upp  (return low if upp <= low)
 */
int getRngInt(int seed = 0, int low = 0, int upp = 9);

/*************************************************************************************************/

/*
 * Generate a pseudo-random boolean value
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param trueFrac: "probability" the the function returns true
 */
bool getRngBool(int seed, double trueFrac = 0.5) { return getRngDouble(seed, 0.0, 1.0) <= trueFrac; }

/*************************************************************************************************/

/*
 * Compute a pseudorandom matrix with elements on a specified range.
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param nRows: number of rows in the output data
 * @param nCols: number of columns in the output data
 * @param[opt] low: lower bound on values in the data  (default: 0.0)
 * @param[opt] upp: upper bound on values in the data  (default: 1.0)
 * @return: a random Eigen array of doubles, using a seed to ensure repeatable results.
 */
Eigen::MatrixXd getRngMatrixXd(int seed, int nRows, int nCols,
                               double low = 0.0, double upp = 1.0);

/*************************************************************************************************/

/*
 * Generate a pseudorandom matrix with elements on a specified range and a specified rank.
 * This is used to test functions that operate on rank-deficient matricies.
 * Note: the values in the matrix are on the order of one, but not strictly bounded.
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param nRows: number of rows in the output data (must be positive)
 * @param nCols: number of columns in the output data (must be positive)
 * @param nRank: the maximum rank that is allowed in the returned matrix
 * @return: a randomly generated matrix that is at most of rank nRank
 */
Eigen::MatrixXd getRngMatrixXdRanked(int seed, int nRows, int nCols, int nRank);

/*************************************************************************************************/

/*
 * Generate a pseudorandom matrix with elements on a specified range and a specified rank.
 * This is used to test functions that operate on rank-deficient matricies.
 * Note: the values in the matrix are on the order of one, but not strictly bounded.
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param nRows: number of elements in the output data (must be positive)
 * @param[opt] low: lower bound on values in the data  (default: 0.0)
 * @param[opt] upp: upper bound on values in the data  (default: 1.0)
 * @return: a randomly generated vector
 */
Eigen::VectorXd getRngVectorXd(int seed, int nRows, double low = 0.0, double upp = 1.0) {
  return getRngMatrixXd(seed, nRows, 1, low, upp).col(0);
}

/*************************************************************************************************/

/*
 * Generate a pseudorandom matrix with elements on a specified range and a specified rank.
 * This is used to test functions that operate on rank-deficient matricies.
 * Note: the values in the matrix are on the order of one, but not strictly bounded.
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then seed is ignored
 * @param low: lower bound on data
 * @param upp: upper bound on data
 * @return: data such that: low <= data <= upp
 */
Eigen::ArrayXd getRngArrBndXd(int seed, const Eigen::ArrayXd& low, const Eigen::ArrayXd& upp);

/*************************************************************************************************/

/*
 * Compute a pseudo-random joint position, drawn from a uniform distribution between the bounds
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then RNG seed is not updated
 * @param qLow: lower bound for each joint
 * @param qUpp: upper bound for each joint
 * @return: joint angles bettween qLow and qUpp.  (return qLow if inconsistent input)
 */
KDL::JntArray getRngBoundedJoints(int seed, const KDL::JntArray& qLow, const KDL::JntArray& qUpp);

/*************************************************************************************************/
/*
 * Compute a pseudo-random joint position, drawn from a uniform distribution between the bounds
 * @param seed: seed to pass to the RNG on each call
 *              if seed == 0, then RNG seed is not updated
 * @param qNom: nominal set of joint angles
 * @param delta: maximum perturbation magnitude along each dimension
 * @param qLow: lower bound for each joint
 * @param qUpp: upper bound for each joint
 * @return: joint angles near qNom and between qLow and qUpp.  (return qNom if inconsistent input)
 */
KDL::JntArray getNearbyJoints(int seed, const KDL::JntArray& qNom, double delta,
                              const KDL::JntArray& qLow, const KDL::JntArray& qUpp);

/*************************************************************************************************/

}  // namespace rng_util
}  // namespace sns_ik

#endif // SNS_IK_RNG_UTILITIES_H_
