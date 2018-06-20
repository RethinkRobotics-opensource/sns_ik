function [ddq, sData, exitCode] = snsIk_acc_QP_mt(ddqLow, ddqUpp, ddxGoalData, JData, dJdqData, alpha)
% [ddq, sData, exitCode] = snsIk_vel_QP(ddqLow, ddqUpp, ddxGoalData, JData, dJdqData, alpha)
%
% This function implements the generic multi-task version of SNS-IK acceleration
% solver using QP.
%
% INPUTS:
%   ddqLow: lower limit for JDataoint accelerations
%   ddqUpp: upper limit for JDataoint accelerations
%   ddxGoalDataData: task-space accelerations of all the tasks in the order of priorities
%   JData: Jacobians of all the task in corresponding order
%   dJdqData: dJ * dq of all the task in corresponding order
%   alpha: trade-off between minimum-joint velocity and maximum task scale
%      default: alpha = 1e-3
%      small alpha: favor solutions that maximize task scale
%      large alpha: favor minimum-joint velocity solutions
%
% OUTPUTS:
%   ddq = Joint acceleration solution with maximum task scale factor
%   sData = task scale factors [0, 1] for all the tasks
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

if nargin < 5
   alpha = 1e-3;
end

% initialize variables
numTask = size(JData,1);
nJnt = numel(ddqLow);
ddqData = cell(numTask,1);
sData = zeros(numTask,1);
resData = cell(numTask,1);
tol = 1e-6;

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

    nDecVar = nJnt + 1;  % [dq; -s]

    % set the objective function
    problem.H = blkdiag(alpha * eye(nJnt), 1.0 / alpha);
    problem.f = zeros(nDecVar, 1);

    % set the constraints:
    problem.lb = [ddqLow; 0];
    problem.ub = [ddqUpp; 1];
    problem.Aineq = [];
    problem.bineq = [];

    % equality constraint for each task (Ai*dq = bi)
    Ai = [Ji, ddxGoali];
    bi = ddxGoali - dJdqi;

    % construct an augmented equality constraints for all the tasks
    % Aa*dq = ba
    if (iTask == 1)
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
    dqiTmp = zSoln(1:nJnt);

    % if the solution becomes too small, the solution is invalid.
    % the current task becomes invalid
    if(norm(dqiTmp) > tol)
        ddqi = ddqiTmp;
        si = 1-zSoln(end);

        % residual
        resi = Ji*ddqi - si*ddxGoali + dJdqi;

        AiBar = [Ji, zeros(nTask,1)];
        biBar = si*ddxGoali - dJdqi + resi;

        % update the previous Aa and ba with residuals
        AaPrev = [AaPrev; AiBar];
        baPrev = [baPrev; biBar];

    else
        % keep the previous solution
        ddqi = ddqData{iTask-1};
        si = 0;
        resi = Ji*ddqi - si*ddxGoali + dJdqi;
        exitCode = 1;
    end

    % store the solutions
    sData(iTask) = si;
    ddqData{iTask} = ddqi;
    resData{iTask} = resi;
end

ddq = ddqi;
end
