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

% This script sets up an runs a set of unit tests on the acceleration IK solvers.

clc; clear;

% shared test parameters;
optTol = 1e-4;
cstTol = 1e-6;
fid = 1;
nTest = 100;

%% SNS-IK acceleration solver, QP implementation
alpha = 1e-4;  % trade-off between minimum joint speed and max task scale
solver = @(dqLow, dqUpp, dx, J, dJdq)( snsIk_acc_QP(dqLow, dqUpp, dx, J, dJdq, alpha) );
result.QP = runTest_snsIk_acc(solver, nTest, optTol, cstTol, fid);

%% SNS-IK solver, Rethink Robotics revised algorithm (Andy Park)
result.rr = runTest_snsIk_acc(@snsIk_acc_rr, nTest, optTol, cstTol, fid);

%% SNS-IK solver, Rethink Robotics revised optimal algorithm (Andy Park)
result.opt = runTest_snsIk_acc(@snsIk_acc_opt, nTest, optTol, cstTol, fid);

%% SNS-IK solver, Rethink Robotics revised algorithm with a secondary configuration space task (Andy Park)
result.rrCS = runTest_snsIk_acc_cs(@snsIk_acc_rr_cs, nTest, optTol, cstTol, fid);
fprintf('average primary task scale factor: %.4f\n', mean(result.rrCS.sData));
fprintf('average secondary task scale factor: %.4f\n', mean(result.rrCS.sCSData));
