function AIC = AIC_RSS(data,data_fit,lambda)
%% Compute AIC from Cluster Assigments using RSS and n_change_points 
% Argha Bandyopadhyay
% Copy of BIC_RSS function, but with penalty term scaled by 2 instead of
% ln(npts)
%
% -------------------------------------------------------------------------
% Overview: 
% ---------
% Compute AIC from residual sum of squares fit where penalty scales with
% number of changePoints in the sequence. Least Absolute Devitations
% Regresssion (L1) is used rather than Least Squares Regression (L2), as it
% is more robust. 
%
% Input Variables:
% ----------------
% data = data to cluster. data is a 1D time series
%
% data_fit = assigment of data_fit  into N labels descibed by cluster means
%
% sigmaNoise = estimated standard deviation of noise from a gaussian
%   distribution. See "estimateNoise.m" for more information. 
%
% Output Variables:
% ----------------
% AIC = AIC value computed for the cluster assignment

%% Compute AIC_RSS
% initalize and check variables
if ~exist('data','var'); disp('Error in AIC_RSS: Need Data to Analyze'); return; end
if ~exist('data_fit','var') || isempty(data_fit); data_fit = ones(length(data),1); end
if length(data) ~= length(data_fit); disp('Error in AIC_RSS: Inputs must be the same length'); return; end
if ~exist('lambda','var') || isempty(lambda)
    lambda = 1; 
end
n_data_points = length(data);
n_change_points = length(find(diff(data_fit) ~= 0));
n_states = length(unique(data_fit));
if data == data_fit
    AIC = (n_change_points + n_states) * 2 * lambda;
else
    AIC = n_data_points * log(sum((data-data_fit).^2)/ n_data_points) +...
        (n_change_points + n_states) * 2 * lambda;
end

