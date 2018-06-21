function [ddq, sData, exitCode] = snsIk_acc_rr_mt(ddqLow, ddqUpp, ddxGoalData, JData, dJdqData)
% [ddq, sData, exitCode] = snsIk_acc_rr_mt(ddqLow, ddqUpp, ddxGoalDataData, JData, dJdqData)
%
% This function implements a generic multi-task version of the SNS-IK
% acceleration solver that tries to meet multiple obJDataective terms
% in a prioritized manner.
%
% INPUTS:
%   ddqLow (nJnt x 1): lower limit for joint accelerations
%   ddqUpp (nJnt x 1): upper limit for joint accelerations
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
%
% This function is a modification of the multi-task SNS-IK algorithm that is
% presented in the original SNS-IK papers. It was developed by Andy Park at
% Rethink Robotics in June 2018.

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

% TODO: input validation

% get the number of tasks
nTask = size(JData,1);
WData = cell(nTask,1);
sData = zeros(nTask,1);
ddqData = cell(nTask,1);
ddqNullData = cell(nTask,1);

% initialization
exitCode = 1;
tol = 1e-6;

nJnt = numel(ddqLow);
I = eye(nJnt);
Pi = I;
ddqi = zeros(nJnt,1);

for iTask = 1:nTask

    % get i-th task Jacobian
    Ji = JData{iTask};
    nddxGoal = size(Ji,1);

    % get i-th task velocity
    ddxGoali = ddxGoalData{iTask};

    % get i-th velocity product
    dJdqi = dJdqData{iTask};

    % update variables for previous projection matrix and solution
    PiPrev = Pi;
    ddqiPrev = ddqi;

    % initialize variables for i-th task
    Wi = eye(nJnt);
    ddqNulli = zeros(nJnt, 1);
    PiBar = PiPrev;
    si = 1.0;
    siStar = 0.0;

    drop_correction_term = false;
    limitExceeded = true;
    cntLoop = 1;
    while limitExceeded == true
        limitExceeded = false;

        % compute a solution without task scale factor
        PiHat = (I - pinv(Ji*PiBar, tol)*Ji)*pinv((I - Wi)*PiPrev, tol);

        if (drop_correction_term)
            ddqi = ddqiPrev + pinv(Ji*PiBar, tol) * (ddxGoali - dJdqi) + ...
                PiHat*(ddqNulli - ddqiPrev);
        else
            ddqi = ddqiPrev + pinv(Ji*PiBar, tol) * (ddxGoali - dJdqi - Ji*ddqiPrev) + ...
                PiHat*(ddqNulli - ddqiPrev);
        end

        % check whether the solution violates the limits
        if any(ddqi < (ddqLow - tol)) || any(ddqi > (ddqUpp + tol))
            limitExceeded = true;
        end

        % compute scale factor
        a = pinv(Ji*PiBar, tol) * ddxGoali;
        b = ddqi - a;

        marginL = ddqLow - b;
        marginU = ddqUpp - b;
        sMax = zeros(nJnt, 1);
        for iJnt=1:nJnt
            if Wi(iJnt,iJnt) == 0
                sMax(iJnt) = inf;
            else
                sMax(iJnt) = FindScaleFactor(marginL(iJnt), marginU(iJnt), a(iJnt));
                % if scale factor is 1 but Joint limit is violated,
                % set a scale factor to a small number so that the
                % corresponding joint will be saturated
                if(sMax(iJnt) == 1 && (ddqi(iJnt) < (ddqLow(iJnt) - tol) || ddqi(iJnt) > (ddqUpp(iJnt) + tol)))
                    sMax(iJnt) = tol;
                end
            end
        end

        [~, JntIdx] = min(sMax);
        taskScale = sMax(JntIdx);

        if (isinf(taskScale))
            % this means that all joints are saturated, so
            % lower priority tasks become infeasible
            taskScale = 0;
        end

        % do the following if the task is feasible
        if (iTask == 1 || taskScale > tol)

            ddqiTmp = taskScale*a + b;
            % try to store maximum scale factor with associated variables
            if (taskScale > siStar || any(ddqiTmp < (ddqLow - tol)) || any(ddqiTmp > (ddqUpp + tol)))
                siStar = taskScale;
                Wistar = Wi;
                ddqNulliStar = ddqNulli;
                PiBarStar = PiBar;
                PiHatStar = PiHat;
            end

            % saturate the most critical joint
            Wi(JntIdx, JntIdx) = 0;
            ddqNulli(JntIdx) = min(max(ddqLow(JntIdx), ddqi(JntIdx)), ddqUpp(JntIdx));

            % update projection matrices
            PiBar = (I - pinv((I - Wi)*PiPrev, tol))*PiPrev;
            PiHat = (I - pinv(Ji*PiBar, tol)*Ji)*pinv((I - Wi)*PiPrev, tol);

            % if rank is below task dimension, terminate the loop and output
            % the current best solution
            if (rank(Ji*PiBar) < nddxGoal)
                si = siStar;
                Wi = Wistar;
                ddqNulli = ddqNulliStar;
                PiBar = PiBarStar;
                PiHat = PiHatStar;

                if (drop_correction_term)
                    ddqi = ddqiPrev + pinv(Ji*PiBar, tol) * (si*ddxGoali - dJdqi) + ...
                        PiHat*(ddqNulli - ddqiPrev);

                    limitExceeded = false;

                    % if the solution for lower priority tasks violates the limits, ignore the solution
                    if (iTask > 1 && (any(ddqi > ddqUpp + tol) || any(ddqi < ddqLow - tol)))
                        si = 0;
                        Wi = zeros(nJnt, nJnt);
                        ddqi = ddqiPrev;
                    end
                else
                    ddqi = ddqiPrev + pinv(Ji*PiBar, tol) * (si*ddxGoali - dJdqi - Ji*ddqiPrev) + ...
                        PiHat*(ddqNulli - ddqiPrev);

                    % if the solution violates the limits, drop the
                    % correction term and run it again
                    if (any(ddqi > ddqUpp + tol) || any(ddqi < ddqLow - tol))
                        drop_correction_term = true;
                        cntLoop = 0;
                        Wi = eye(nJnt);
                        ddqNulli = zeros(nJnt, 1);
                        PiBar = PiPrev;
                        si = 1.0;
                        siStar = 0.0;
                    else
                        limitExceeded = false;
                    end
                end
            end

        else % if the current task is infeasible
            si = 0;
            Wi = zeros(nJnt, nJnt);
            ddqi = ddqiPrev;
            limitExceeded = false;
        end

        cntLoop = cntLoop + 1;
    end

    % only if the current task is feasible, update nullspace projection
    if (si > 0)
        Pi = PiPrev - pinv(Ji*PiPrev, tol)*(Ji*PiPrev);
    end

    % store data
    WData{iTask} = Wi;
    sData(iTask) = si;
    ddqData{iTask} = ddqi;
    ddqNullData{iTask} = ddqNulli;
end

%-- end of algorithm
ddq = ddqi;

end
