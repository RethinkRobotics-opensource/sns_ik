function result = runTest_snsIk_acc(solver, nTest, optTol, cstTol, fid)
% result = runTest_snsIk_acc(solver, nTest, optTol, cstTol, fid)
%
% This function runs a unit test on a candidate acceleration IK solver
%
% INPUTS:
%   solver: acceleration Ik solver
%       [ddq, s, exitCode] = solver(ddqLow, ddqUpp, ddx, J)
%       IN: ddqLow = lower bound on joint velocity
%       IN: ddqUpp = upper bound on joint velocity
%       IN: ddx = task velocity (goal)
%       IN: J = task jacobian
%       IN: dJddq = task jacobian rate * joint angle rate
%       OUT: ddq = joint velocity (solution)
%       OUT: s = task scale
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
%

% run the tests
nPass = 0;
nFail = 0;
solveTime = zeros(nTest, 1);
for iTest = 1:nTest
    fprintf(fid, '--------------- TEST %d ---------------', iTest);

    % initialize the RNG seed for repeatable tests
    rng(iTest + 4283743);

    % set up the problem dimensions
    nJnt = randi(8);
    nTask = randi(nJnt);

    % set up the problem itself
    J = randn(nTask, nJnt);  % jacobian
    dJdq = randn(nTask, 1);  % jacobian rate * joint velocity
    ddqLow = -0.1 - rand(nJnt, 1);
    ddqUpp = 0.1 + rand(nJnt, 1);

    % compute the unscaled solution
    ddqTmp = ddqLow + (ddqUpp - ddqLow) .* rand(nJnt, 1);
    for iJnt = 1:nJnt
        if rand(1) < 0.2
            if rand(1) < 0.5
                ddqTmp(iJnt) = ddqLow(iJnt);
            else
                ddqTmp(iJnt) = ddqUpp(iJnt);
            end
        end
    end

    % compute the task velocity, then scale
    scale = 1.0;
    if rand(1) < 0.4;
        scale = 0.1 + 0.8 * rand(1);
    end
    ddxGoal = (J * ddqTmp + dJdq)/ scale;

    % solve the problem:
    testPass = true;
    startTime = tic;
    [ddq, s, exitCode] = solver(ddqLow, ddqUpp, ddxGoal, J, dJdq);
    solveTime(iTest) = toc(startTime);

    % check the solution
    if exitCode ~= 1
        fprintf(fid, '\n  solver failed with exit code %d', exitCode);
        testPass = false;
    end
    if s < scale - optTol
        fprintf(fid, '\n  task scale is sub-optimal!  s (%.6f) < scale (%.6f)', s, scale);
        testPass = false;
    end
    if any(ddq > ddqUpp + cstTol) || any(ddq < ddqLow - cstTol)
        fprintf(fid, '\n  joint limits are violated!  [ddqLow, ddq, ddqUpp]');
        testPass = false;
    end
    taskError = s*ddxGoal - J*ddq - dJdq;
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

end
