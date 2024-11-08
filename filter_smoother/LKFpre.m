function [X_ests,dx_ests,P_ests,P_preds,bs] = LKFpre(X_est0,dx_est0,P_0,Hs,Qs,Rs,zs,X_preds,Phis)
%LKFpre takes in an initial state, state deviation, and covariance; 
% observation, process noise, and measurement noise matrices; and 
% precomputed lists of observations, reference trajectory and error 
% state transition matrices in order to perform a Linearized Kalman Filter 
% pass on the precomputed data and initial conditions.
% INPUTS:
%  X_est0 - [m,1] Initial State Estimate
% dx_est0 - [m,1] Initial State Uncertainty Estimate
%     P_0 - [m,m] Initial State Covariance Matrix
%      Hs - [k,m,(1 or n)] Observation Matrix
%      Qs - [m,m,(1 or n)] Process Noise Matrix
%      Rs - [m,m,(1 or n)] Measurement Noise Matrix
%      zs - [n,k] Precomputed Observations
% X_preds - [n,m] Precomputed Reference States
%    Phis - [m,m,(1 or n)] Precomputed (Error-)State Transition Matrices
% OUTPUT:
%  X_ests - [n,m] Filtered State Estimates
% dx_ests - [n,m] Filtered State Uncertainty Estimates
%  P_ests - [m,m,n] Estimated Covariance Matrices
% P_preds - [m,m,n] Predicted Covariance Matrices
%      bs - [n,k] Filter Innovations (Predicted Deviations from Observations)

[N,K] = size(zs);
M = length(X_est0);

% Initialize filter
X_est = X_est0;
dx_est = dx_est0;
P_est = P_0;

% Preallocate output
X_ests = zeros(N,M);
dx_ests = zeros(N,M);
P_ests = zeros(M,M,N);
P_preds = zeros(M,M,N);
bs = zeros(N,K);

% Set initial output
X_ests(1,:) = X_est;
dx_ests(1,:) = dx_est0;
P_ests(:,:,1) = P_0;
P_preds(:,:,1) = P_0;

% Handle Single or Multidimensional Cases
% Precompute logic to see if the matrices change
vH = ~ismatrix(Hs);
vQ = ~ismatrix(Qs);
vR = ~ismatrix(Rs);
vPhi = ~ismatrix(Phis);

% Iterate over the remainder of observations
for i = 2:N
    % We set these here to avoid extra indexing later
    if vH % H varies
        H = Hs(:,:,i);
    else  % H constant
        H = Hs;
    end
    if vQ % Q varies
        Q = Qs(:,:,i);
    else  % Q constant
        Q = Qs;
    end
    if vR % R varies
        R = Rs(:,:,i);
    else  % R constant
        R = Rs;
    end
    if vPhi % Phi varies
        Phi = Phis(:,:,i);
    else    % Phi constant (linear system)
        Phi = Phis;
    end

    % Load precomputed obervations and reference
    z = zs(i,:).';
    X_pred = X_preds(i,:).';

    % Time Update
    dx_pred = Phi*dx_est;
    P_pred = Phi*P_est*Phi' + Q;

    % Compute observation deviation
    if allfinite(z)
        b = z - H*X_pred;
    else
        b = H*X_pred;
        wid = "Value:ObservationNaNorInf";
        msg = "NaN or Inf observation detected.";
        warning(wid,msg)
    end

    %Compute Kalman gain matrix
    K = P_pred*H' / (H*P_pred*H' + R);

    % Measurement Update
    dx_est = dx_pred + K*(b-H*dx_pred);
    %P_est = P_pred - K*H*P_pred;
    X_est = X_pred + dx_est;

    % Implement Joseph Stabilized Formulation for Covariance Update
    J = speye(M) - K*H;
    if allfinite(R)
        P_est = J*P_pred*J.' + K*R*K.';
    else
        % The measurement noise is infinite or not available
        % Thus K -> 0 and K*R*K' -> 0
        P_est = J*P_pred*J.';
        wid = "Value:MeasurementNoiseNaNorInf";
        msg = "NaN or Inf measurement noise detected.";
        warning(wid,msg)
    end

    % Ensure P_est is symmetric
    P_est = (P_est.' + P_est)/2;

    % Save outputs into arrays
    X_ests(i,:) = X_est.';
    dx_ests(i,:) = dx_est.';
    P_ests(:,:,i) = P_est;
    P_preds(:,:,i) = P_pred;
    bs(i,:) = b.';
end

end