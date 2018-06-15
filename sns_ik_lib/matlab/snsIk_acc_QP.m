function [ddq, s, exitCode] = snsIk_acc_QP(ddqLow, ddqUpp, ddxGoal, J, dJdq, alpha)
% [ddq, s, exitCode] = snsIk_acc_QP(ddqLow, ddqUpp, ddxGoal, J, dJdq, alpha)
%
% This function implements the QP version of the SNS-IK acceleration solver
%
% INPUTS:
%   ddqLow: lower limit for joint acceleration
%   ddqUpp: upper limit for joint acceleration
%   ddxGoal: task-space acceleration
%   J: jacoabian mapping joint space to task space
%   dJdq = dJ * dq
%       dJ = time-derivative of the task jacobian
%       dq = current joint velocity
%   alpha: trade-off between minimum-joint acceleration and maximum task scale
%      default: alpha = 1e-3
%      small alpha: favor solutions that maximize task scale
%      large alpha: favor minimum-joint acceleration solutions
%
% OUTPUTS:
%   ddq = joint acceleration solution with maximum task scale factor
%   s = task scale factor [0, 1]
%
% NOTES:
%   s * ddxGoal = J * ddq + dJ * dq;
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

% TODO: input validation
% TODO: return optimization status

if nargin < 6
   alpha = 1e-3;
end

% problem dimensions
[~, nJnt] = size(J);
nDecVar = nJnt + 1;  % [ddq; -s]

% set the objective function
problem.H = blkdiag(alpha * eye(nJnt), 1.0 / alpha);
problem.f = zeros(nDecVar, 1);

% set the constraints:
problem.lb = [ddqLow; 0];
problem.ub = [ddqUpp; 1];
problem.Aeq = [J, ddxGoal];
problem.beq = ddxGoal - dJdq;
problem.Aineq = [];
problem.bineq = [];

% solve the problem
problem.options = optimset('Display','off');
problem.solver = 'quadprog';
[zSoln, ~, exitCode] = quadprog(problem);
ddq = zSoln(1:nJnt);
s = 1 - zSoln(end);

end
