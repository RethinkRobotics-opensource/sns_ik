function [dq, s, exitCode] = snsIk_vel_QP(dqLow, dqUpp, dxGoal, J, alpha)
% [dq, s, exitCode] = snsIk_vel_QP(dqLow, dqUpp, dxGoal, J, alpha)
%
% This function implements the basic QP version of the SNS-IK velocity
% solver.
%
% INPUTS:
%   dqLow (nJnt x 1): lower limit for joint velocities
%   dqUpp (nJnt x 1): upper limit for joint velocities
%   dxGoal (nJnt x 1): task-space velocity
%   J (ndxGoal x nJnt): jacoabian mapping joint space to task space
%   alpha: trade-off between minimum-joint velocity and maximum task scale
%      default: alpha = 1e-3
%      small alpha: favor solutions that maximize task scale
%      large alpha: favor minimum-joint velocity solutions
%
% OUTPUTS:
%   dq (nJnt x 1): joint velocity solution with maximum task scale factor
%   s = task scale factor [0, 1]
%
% NOTES:
%   s * dxGoal = J * dq;

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

% problem dimensions
[nTask, nJnt] = size(J);
nDecVar = nJnt + 1;  % [dq; -s]

% set the objective function
problem.H = blkdiag(alpha * eye(nJnt), 1.0 / alpha);
problem.f = zeros(nDecVar, 1);

% set the constraints:
problem.lb = [dqLow; 0];
problem.ub = [dqUpp; 1];
problem.Aeq = [J, dxGoal];
problem.beq = dxGoal;
problem.Aineq = [];
problem.bineq = [];

% solve the problem
problem.options = optimset('Display','off');
problem.solver = 'quadprog';
[zSoln, ~, exitCode] = quadprog(problem);
dq = zSoln(1:nJnt);
s = 1-zSoln(end);

end
