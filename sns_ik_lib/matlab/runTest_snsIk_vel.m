function result = runTest_snsIk_vel(solver, nTest, optTol, cstTol, fid)
% result = runTest_snsIk_vel(solver, nTest, optTol, cstTol, fid)
%
% This function runs a unit test on a candidate velocity IK solver
%
% INPUTS:
%   solver: velocity Ik solver
%       [dq, s, exitCode] = solver(dqLow, dqUpp, dx, J)
%       IN: dqLow = lower bound on joint velocity
%       IN: dqUpp = upper bound on joint velocity
%       IN: dx = task velocity (goal)
%       IN: J = task jacobian
%       OUT: dq = joint velocity (solution)
%       OUT: s = task scale
%       OUT: exitCode  (1 == success)
%   nTest: how many tests to run?
%   optTol: optimality tolerance
%   cstTol: feasibility tolerance
%   fid: output stream for test logs (1 == command prompt)
%
% OUTPUTS:
%   result: struct with test results
%

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
    rng(iTest + 9872342);

    % set up the problem dimensions
    nJnt = randi(8);
    nTask = randi(nJnt);

    % set up the problem itself
    J = randn(nTask, nJnt);  % jacobian
    dqLow = -0.1 - rand(nJnt, 1);
    dqUpp = 0.1 + rand(nJnt, 1);

    % compute the unscaled solution
    dqTmp = dqLow + (dqUpp - dqLow) .* rand(nJnt, 1);
    for iJnt = 1:nJnt
        if rand(1) < 0.2
            if rand(1) < 0.5
                dqTmp(iJnt) = dqLow(iJnt);
            else
                dqTmp(iJnt) = dqUpp(iJnt);
            end
        end
    end

    % compute the task velocity, then scale
    scale = 1.0;
    if rand(1) < 0.4;
        scale = 0.1 + 0.8 * rand(1);
    end
    dxGoal = J * dqTmp / scale;

    % solve the problem:
    testPass = true;
    startTime = tic;
    [dq, s, exitCode] = solver(dqLow, dqUpp, dxGoal, J);
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
    if any(dq > dqUpp + cstTol) || any(dq < dqLow - cstTol)
        fprintf(fid, '\n  joint limits are violated!  [dqLow, dq, dqUpp]');
        testPass = false;
    end
    taskError = s*dxGoal - J*dq;
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
