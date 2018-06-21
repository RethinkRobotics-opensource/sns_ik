function result = runTest_snsIk_vel_mt(solver, nTest, optTol, cstTol, fid)
% result = runTest_snsIk_vel_mt(solver, nTest, optTol, cstTol, fid)
%
% This function runs a unit test on a candidate velocity IK solver
% that supports multiple task objective terms.
%
% INPUTS:
%   solver: velocity Ik solver
%       [dq, sData, exitCode] = solver(dqLow, dqUpp, dxData, JData)
%       IN: dqLow (nJnt x 1) = lower bound on joint velocity
%       IN: dqUpp (nJnt x 1) = upper bound on joint velocity
%       IN: dxGoalData (cell array) = task velocities
%       IN: JData (cell array) = task jacobians
%       OUT: dq (nJnt x 1) = joint velocity (solution)
%       OUT: sData (nTask x 1) = task scales
%       OUT: exitCode  (1 == success)
%   nTest: how many tests to run?
%   optTol: optimality tolerance
%   cstTol: feasibility tolerance
%   fid: output stream for test logs (1 == command prompt)
%
% OUTPUTS:
%   result: struct with test results

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

% run the tests
nPass = 0;
nFail = 0;
solveTime = zeros(nTest, 1);
sDataTotal = cell(nTest,1);

for iTest = 1:nTest
    fprintf(fid, '--------------- TEST %d ---------------', iTest);

    % initialize the RNG seed for repeatable tests
    rng(iTest + 5465454);

    % set up the problem dimensions
    nJnt = 2 + randi(8);

    dqLow = -0.1 - rand(nJnt, 1);
    dqUpp = 0.1 + rand(nJnt, 1);

    nTask = 2 + randi(3); % number of tasks
    ndxGoal = cell(nTask,1);
    scale = ones(nTask,1);
    JData = cell(nTask,1);
    dqTmp = zeros(nJnt,nTask);
    dxGoalData = cell(nTask,1);

    for iTask = 1:nTask
        ndxGoal{iTask} = randi(nJnt-2); % task space dimension

        % set up the problem itself
        JData{iTask} = randn(ndxGoal{iTask}, nJnt);  % jacobian

        % compute the unscaled solution
        dqTmp(:,iTask) = dqLow + (dqUpp - dqLow) .* rand(nJnt, 1);

        for iJnt = 1:nJnt
            if rand(1) < 0.2
                if rand(1) < 0.5
                    dqTmp(iJnt,iTask) = dqLow(iJnt);
                else
                    dqTmp(iJnt,iTask) = dqUpp(iJnt);
                end
            end
        end

        % compute the task velocity, then scale
        if rand(1) < 0.4
            scale(iTask) = 0.1 + 0.8 * rand(1);
        end
        dxGoalData{iTask} = JData{iTask} * dqTmp(:,iTask) / scale(iTask);

    end

    % solve the problem:
    testPass = true;
    startTime = tic;
    [dq, sData, exitCode] = solver(dqLow, dqUpp, dxGoalData, JData);
    solveTime(iTest) = toc(startTime);
    sDataTotal{iTest} = sData;
    sPrimaryTask = sDataTotal{iTest}(1);
    scalePrimaryTask = scale(1);

    % check the solution
    if exitCode ~= 1
        fprintf(fid, '\n  solver failed with exit code %d', exitCode);
        testPass = false;
    end
    if sPrimaryTask < scalePrimaryTask - optTol
        fprintf(fid, '\n  task scale is sub-optimal!  s (%.6f) < scale (%.6f)', sPrimaryTask, scalePrimaryTask);
        testPass = false;
    end
    if any(dq > dqUpp + cstTol) || any(dq < dqLow - cstTol)
        fprintf(fid, '\n  joint limits are violated!  [dqLow, dq, dqUpp]');
        testPass = false;
    end
    taskError = sPrimaryTask*dxGoalData{1} - JData{1}*dq;
    if any(abs(taskError) > cstTol)
        fprintf(fid, '\n  task error (%.6f) exceeds tolerance!', max(abs(taskError)));
        testPass = false;
    end
    if testPass
        nPass = nPass + 1;
        fprintf(fid, '  --> pass\n');
    else
        nFail = nFail + 1;
        fprintf(fid, '\n  --> fail\n');
    end
end
fprintf(fid, '\n');
fprintf(fid, '------------------------------------\n');
fprintf(fid, 'nPass: %d  --  nFail: %d\n', nPass, nFail);
fprintf(fid, '------------------------------------\n');

result.nPass = nPass;
result.nFail = nFail;
result.solveTime = solveTime;
result.sDataTotal = sDataTotal;
end
