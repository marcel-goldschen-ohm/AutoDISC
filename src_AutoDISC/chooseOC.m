function optimal_OC = chooseOC(SNR, numsamples, slope, intercept)
%%  Choose OC
%   Author: Argha Bandyopadhyay
%   Contact: argha.bandyopadhyay@utexas.edu
%
%   Published in:
%   Bandyopadhyay and Goldschen-Ohm, 2021
%   -----------------------------------
%   **Overview:**
%   Chooses the optimal OC (either BIC_GMM or BIC_RSS) based
%   on the SNR and number of samples in the time series.
%   Choice mediated by two-dimensional decision boundary
%   (either default or user-specified) - all points above
%   the decision boundary are GMM-optimal, all points below
%   the decision boundary are RSS-optimal.
%
%   -----------------------------------
%   **I/O:**
%   Inputs:
%       1. SNR: double which represents state separation-to-noise 
%       ratio of input trace
%       2. numsamples: integer number of points (length) of input
%       trace
%       3. db_input: struct with slope and intercept of
%       decision boundary
%           db_input.slope: double which represents decision
%           boundary slope - default = -0.74
%           db_input.intercept: double which represents
%           decision boundary intercept - default = 5.8
%   Outputs:
%       1. optimal_OC: string of either "BIC_GMM" or
%       "BIC_RSS" depending on trace's relation to decision
%       boundary
%   -----------------------------------



%   Check variables and set defaults
if ~exist('SNR', 'var') || ~isnumeric(SNR)
    disp("Error: SNR of time series not provided in correct format (numeric).");
end

if ~exist('numsamples', 'var') || ~isnumeric(numsamples)
    disp("Error: Length of time series not provided in correct format (numeric).");
end

if ~exist('slope', 'var') || ~isnumeric(slope)
    disp("Error: Slope of decision boundary not provided in correct format (numeric).");
end

if ~exist('intercept', 'var') || ~isnumeric(intercept)
    disp("Error: Intercept of decision boundary not provided in correct format (numeric).");
end


%   Decision boundary relates numsamples to SNR in the
%   following formula:
%   log10(numsamples) = db_input.slope*SNR + db_input.intercept
decision_boundary = slope*SNR + intercept;

%   if log10(numsamples) above boundary, then BIC_GMM is optimal
%   if log10(numsamples) below or on boundary, then BIC_RSS is optimal
if log10(numsamples) > decision_boundary
    optimal_OC = "AIC_GMM";
else
    optimal_OC = "BIC_RSS";
end
end
