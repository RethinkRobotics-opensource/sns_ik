function result = runTest_snsIk_acc_mt(solver, nTest, optTol, cstTol, fid)
% result = runTest_snsIk_acc_mt(solver, nTest, optTol, cstTol, fid)
%
% This function runs a unit test on a candidate acceleration IK solver
% that supports multiple task objective terms.
%
% INPUTS:
%   solver: acceleration Ik solver
%       [ddq, sData, exitCode] = solver(ddqLow, ddqUpp, ddxData, JData)
%       IN: ddqLow (nJnt x 1)= lower bound on joint acceleration
%       IN: ddqUpp (nJnt x 1)= upper bound on joint acceleration
%       IN: ddxGoalData (cell array)= task accelerations
%       IN: JData (cell array)= task jacobians
%       IN: dJdqData (cell array)= the product of task jacobian derivative and joint accelerations
%       OUT: ddq (nJnt x 1)= joint acceleration (solution)
%       OUT: sData (nTask x 1)= task scales
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
    rng(iTest + 4283743);

    % set up the problem dimensions
    nJnt = 2 + randi(7);
    ddqLow = -0.1 - rand(nJnt, 1);
    ddqUpp = 0.1 + rand(nJnt, 1);

    nTask = 2 + randi(4); % number of tasks
    nddxGoal = cell(nTask,1);
    scale = ones(nTask,1);
    JData = cell(nTask,1);
    dJdqData = cell(nTask,1);
    ddqTmp = zeros(nJnt,nTask);
    ddxGoal = cell(nTask,1);

    for iTask = 1:nTask
        nddxGoal{iTask} = randi(nJnt); % task space dimension

        % set up the problem itself
        JData{iTask} = randn(nddxGoal{iTask}, nJnt);  % jacobian
        dJdqData{iTask} = randn(nddxGoal{iTask}, 1);  % jacobian rate * joint velocity

        % compute the unscaled solution
        ddqTmp(:,iTask) = ddqLow + (ddqUpp - ddqLow) .* rand(nJnt, 1);
        for iJnt = 1:nJnt
            if rand(1) < 0.2
                if rand(1) < 0.5
                    ddqTmp(iJnt,iTask) = ddqLow(iJnt);
                else
                    ddqTmp(iJnt,iTask) = ddqUpp(iJnt);
                end
            end
        end

        % compute the task velocity, then scale
        if rand(1) < 0.4
            scale(iTask) = 0.1 + 0.8 * rand(1);
        end
        ddxGoal{iTask} = (JData{iTask} * ddqTmp(:,iTask) + dJdqData{iTask})/ scale(iTask);

    end

    % solve the problem:
    testPass = true;
    startTime = tic;
    [ddq, sData, exitCode] = solver(ddqLow, ddqUpp, ddxGoal, JData, dJdqData);
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
    if any(ddq > ddqUpp + cstTol) || any(ddq < ddqLow - cstTol)
        fprintf(fid, '\n  joint limits are violated!  [ddqLow, ddq, ddqUpp]');
        testPass = false;
    end
    taskError = sPrimaryTask*ddxGoal{1} - JData{1}*ddq - dJdqData{1};
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
