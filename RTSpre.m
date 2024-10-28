function [X_sms,dx_sms,P_sms] = RTSpre(X_preds,dx_ests,Phis,P_ests,P_preds)
%RTSpre takes the precalculated reference trajectory and filtered estimated
% uncertainties along with the error state transition matrices, covariance
% estimates and predictions from the Kalman Filter pass. Performs a single
% pass of a Rauch–Tung–Striebel Smoother. 
% INPUTS:
%  X_preds - [n,m] Precomputed Reference States
%  dx_ests - [n,m] Filtered State Uncertainty Estimates
%     Phis - [m,m,n] Precomputed (Error-)State Transition Matrices
%   P_ests - [m,m,n] Estimated Covariance Matrices
%  P_preds - [m,m,n] Predicted Covariance Matrices
% OUTPUT:
%      X_sm - [1x1] Smoothed States
%     dx_sm - [1x1] Smoothed State Uncertainty/Errors
%      P_sm - [m,m,n] Smoothed Covariance Matrices

[N,M] = size(X_preds);

% Preallocate Outputs
dx_sms = zeros(N,M);
P_sms = zeros(M,M,N);

% Initialize Smoother
dx_sm = dx_ests(N,:).';
P_sm = P_ests(:,:,N);

% Save last iteration to outputs
dx_sms(N,:) = dx_sm.';
P_sms(:,:,N) = P_sm;

for i = N-1:-1:1
    % Read matrices
    P_1 = P_ests(:,:,i);
    P_2 = P_preds(:,:,i+1);
    Phi = Phis(:,:,i+1);
    dx_est = dx_ests(i,:).';

    % Perform Smoothing
    S = P_1*Phi.' / P_2;
    dx_sm = dx_est + S*(dx_sm - Phi*dx_est);
    P_sm = P_1 + S*(P_sm-P_2)*S.';

    % Save outputs into arrays
    dx_sms(i,:) = dx_sm.';
    P_sms(:,:,i) = P_sm;
end

X_sms = X_preds + dx_sms;

end