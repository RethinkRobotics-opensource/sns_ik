function [dq, s, exitCode] = snsIk_vel_basic(dqLow, dqUpp, dx, J)
% [dq, s, exitCode] = snsIk_vel_basic(dqLow, dqUpp, dx, J)
%
% This function implements the basic version of the SNS-IK velocity
% solver.
%
% INPUTS:
%   dqLow: lower limit for joint velocities
%   dqUpp: upper limit for joint velocities
%   dx: task-space velocity
%   J: jacoabian mapping joint space to task space
%
% OUTPUTS:
%   dq = joint velocity solution with maximum task scale factor
%   s = task scale factor [0, 1]
%   exitCode    (1 == success)
%
%
% NOTES:
%   s * dxGoal = J * dq;
%
%  --> This implementation is as close to the standard SNS-IK algorithm as
%      is possible, outlined in the main SNS-IK paper.
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

% This implementation is not reliable, and in many cases will fail to
% compute the correct solution. It is unclear whether this is a bug in
% this implementation, or a bug in the algorithm that is described in the
% paper. See snsIk_vel_rr() for a modified version of this algorithm that
% does pass all of the unit tests.
warning('This implementation does not work!');

[nTask, nJnt] = size(J);

W = eye(nJnt);
dqNull = zeros(nJnt, 1);
s = 1.0;
sStar = 0.0;
maxIter = 5 * nJnt;  % over-estimate of how many iterations to expect
for iter = 1:maxIter
    limitExceeded = false;
    dq = dqNull + pinv(J*W) * (dx - J*dqNull);
    if any(dq < dqLow) || any(dq > dqUpp)
       limitExceeded = true;
    end
    a = pinv(J*W) * dx;
    b = dq - a;
    [taskScale, jntIdx] = getTaskScalingFactor(a, b, dqLow, dqUpp);
    if taskScale > sStar
        sStar = taskScale;
        Wstar = W;
        dqNullStar = dqNull;
    end
    W(jntIdx, jntIdx) = 0;
    dq(jntIdx) = min(max(dqLow(jntIdx), dq(jntIdx)), dqUpp(jntIdx));
    if rank(J*W) < nTask
       s = sStar;
       W = Wstar;
       dqNull = dqNullStar;
       dq = dqNull + pinv(J*W)*(s * dx - J * dqNull);
%        limitExceeded = false;  % unreachable
       break;
    end
    if limitExceeded
        break;
    end
end

% TODO: better exit code
if iter < maxIter
    exitCode = 1;
else
    exitCode = -1;
end

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [taskScale, jntIdx] = getTaskScalingFactor(a, b, dqLow, dqUpp)

% TODO: documentation

sMin = (dqLow - b) ./ a;
sMax = (dqUpp - b) ./ a;
for i=1:length(sMin)
   if sMin(i) > sMax(i)
      tmp = sMin(i);
      sMin(i) = sMax(i);
      sMax(i) = tmp;
   end
end

[sMaxVal, sMaxIdx] = min(sMax);
sMinVal = max(sMin);
if sMinVal > sMaxVal || sMaxVal < 0.0 || sMinVal > 1.0
    taskScale = 0.0;
    warning('Infeasible solution!');
else
    taskScale = sMaxVal;
    jntIdx = sMaxIdx;
end

end
