% Copyright 2018 Rethink Robotics
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Run a simple unit test on the matlab implementation for the various
% velocity IK solvers that are related to the original SNS-IK papers.

clc; clear;

% shared test parameters;
optTol = 1e-4;
cstTol = 1e-6;
fid = 1;
nTest = 100;

%% SNS-IK velocity solver, QP implementation
alpha = 1e-4;  % trade-off between minimum joint speed and max task scale
solver = @(dqLow, dqUpp, dx, J)( snsIk_vel_QP(dqLow, dqUpp, dx, J, alpha) );
result.QP = runTest_snsIk_vel(solver, nTest, optTol, cstTol, fid);

%% SNS-IK solver, Rethink Robotics revised algorithm (Andy Park)
result.rr = runTest_snsIk_vel(@snsIk_vel_rr, nTest, optTol, cstTol, fid);

%% SNS-IK solver, Rethink Robotics revised optimal algorithm (Andy Park)
result.opt = runTest_snsIk_vel(@snsIk_vel_opt, nTest, optTol, cstTol, fid);

%% SNS-IK solver, original SNS-IK algorithm
% % This solver currently does not pass the test
% result.basic = runTest_snsIk_vel(@snsIk_vel_basic, nTest, optTol, cstTol, fid);
