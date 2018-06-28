/** @file test_utilities.hpp
 *
 * @brief a collection of functions for testing Eigen type vector data
 *
 * @author Matthew Kelly
 * @author Andy Park
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

#ifndef SNS_IK_LIB_TEST_UTILITIES_H
#define SNS_IK_LIB_TEST_UTILITIES_H

#include <Eigen/Dense>

namespace sns_ik {
namespace test_util {

/**
 * Check if the elements in vector A and B are equal
 */
void checkEqualVector(const Eigen::VectorXd& A, const Eigen::VectorXd& B, double tol);

/**
 * Check if the elements in a vector are within lower and upper limits
 */
void checkVectorLimits(const Eigen::VectorXd& low, const Eigen::VectorXd& val, const Eigen::VectorXd& upp, double tol);

}  // namespace test_util
}  // namespace sns_ik

#endif //SNS_IK_LIB_TEST_UTILITIES_H