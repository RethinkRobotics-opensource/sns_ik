function [ddq, sData, exitCode] = snsIk_acc_QP_mt(ddqLow, ddqUpp, ddxGoalData, JData, dJdqData)
% [ddq, sData, exitCode] = snsIk_vel_QP(ddqLow, ddqUpp, ddxGoalData, JData, dJdqData)
%
% This function implements the generic multi-task version of SNS-IK acceleration
% solver using QP.
%
% INPUTS:
%   ddqLow (nJnt x 1): lower limit for JDataoint accelerations
%   ddqUpp (nJnt x 1): upper limit for JDataoint accelerations
%   ddxGoalData (cell array): task-space accelerations of all the tasks in the order of priorities
%   JData (cell array): Jacobians of all the task in corresponding order
%   dJdqData (cell array): dJ * dq of all the task in corresponding order
%
% OUTPUTS:
%   ddq (nJnt x 1): Joint acceleration solution with maximum task scale factor
%   sData (nTask x 1): task scale factors [0, 1] for all the tasks
%   exitCode    (1 == success)
%
% NOTES:

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

% TODO: input validation
% TODO: return optimization status
% TODO: see if we can solve all tasks in a single QP problem

%   alpha: trade-off between minimum-joint velocity and maximum task scale
%    default: alpha = 1e-3
%    small alpha: favor solutions that maximize task scale
%    large alpha: favor minimum-joint velocity solutions
alpha = 1e-3;

% initialize variables
numTask = size(JData,1);
nJnt = numel(ddqLow);
ddqData = cell(numTask,1);
sData = zeros(numTask,1);
resData = cell(numTask,1);

% initialize problem
problem.options = optimset('Display','off');
problem.solver = 'quadprog';

for iTask = 1:numTask

    % get i-th task jacobian
    Ji = JData{iTask};
    nTask = size(Ji,1);

    % get i-th task velocity
    ddxGoali = ddxGoalData{iTask};

    % get i-th velocity product
    dJdqi = dJdqData{iTask};

    % equality constraint for each task (Ai*dq = bi)
    Ai = [Ji, ddxGoali];
    bi = ddxGoali - dJdqi;

    % construct an augmented equality constraints for all the tasks
    % Aa*dq = ba
    if (iTask == 1)
        % set problem dimension
        nDecVar = nJnt + 1;  % [dq; -s]

        % set the objective function
        problem.H = blkdiag(alpha * eye(nJnt), 1.0 / alpha);
        problem.f = zeros(nDecVar, 1);

        % set the constraints:
        problem.lb = [ddqLow; 0];
        problem.ub = [ddqUpp; 1];
        problem.Aineq = [];
        problem.bineq = [];

        problem.Aeq = Ai;
        problem.beq = bi;
        AaPrev = [];
        baPrev = [];
    else
        problem.Aeq = [AaPrev; Ai];
        problem.beq = [baPrev; bi];
    end

    % solve the problem
    [zSoln, ~, exitCode] = quadprog(problem);

    % if the solution becomes too small, the solution is invalid.
    % the current task becomes invalid
    if(exitCode == 1)
        ddqi = zSoln(1:nJnt);
        si = 1-zSoln(end);

        % residual
        resi = Ji*ddqi - si*ddxGoali + dJdqi;

        AiBar = [Ji, zeros(nTask,1)];
        biBar = si*ddxGoali - dJdqi + resi;

        % update the previous Aa and ba with residuals
        AaPrev = [AaPrev; AiBar];
        baPrev = [baPrev; biBar];

   elseif (exitCode < 1 && iTask > 1)
        % if the solution for a lower priority task does not exist,
        % the current task becomes infeasible (i.e., si = 0)
        % we will keep the previous solution

        ddqi = ddqData{iTask-1};
        si = 0;
        resi = Ji*ddqi - si*ddxGoali + dJdqi;
        exitCode = 1;

    else
        dqi = zeros(nJnt,1);
        si = 0;
        resi = Ji*dqi - si*dxGoali;
    end

    % store the solutions
    sData(iTask) = si;
    ddqData{iTask} = ddqi;
    resData{iTask} = resi;
end

ddq = ddqi;
end
