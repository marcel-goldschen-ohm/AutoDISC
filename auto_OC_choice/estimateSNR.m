function estimated_SNR = estimateSNR(timeseries, idealization)
%%  Estimate SNR
%   Author: Argha Bandyopadhyay
%   Contact: argha.bandyopadhyay@utexas.edu
%
%   Published in:
%   Bandyopadhyay and Goldschen-Ohm, 2021
%   -----------------------------------
%   **Overview:**
%   Estimates SNR of a time series given a rough
%   idealization of the trace.
%   
%   Workflow:
%       1. Signal Estimation: Binary k-Means clustering on
%       the amplitudes of state transitions found by the
%       idealization. Assumes that low amplitude cluster is
%       representative of overfitting to noise, so high
%       amplitude cluster is representative of true "signal"
%       which is average state separation
%       2. Noise Estimation: Residuals between the
%       idealization and the noisy time series are
%       calculated, with their standard deviation
%       representing the average noise level of the trace
%       3. SNR Calculation: Division of estimated signal by
%       estimated noise
%   -----------------------------------
%   **Requirements:**
%   MATLAB Statistics and Machine Learning Toolbox
%   -----------------------------------
%   **I/O:**
%   Inputs:
%       1. timeseries: vector of noisy data with size
%       [N,1] where N >= 5
%       2. idealization: vector of ideal state amplitude
%       sequence with size [N,1]
%   Outputs:
%       1. estimated_SNR = SNR estimate double
%   -----------------------------------



%   Check variables
if ~exist('timeseries', 'var') || length(timeseries) < 5
    disp("Error: Noisy Time Series not found or too short to analyze.");
end

if ~exist('idealization', 'var')
    disp("Error: Idealization not found");
end

if length(idealization) ~= length(timeseries)
    disp("Error: Length of idealization does not match length of noisy time series");
end

% Do not want to get warnings about kmeans failing to
% converge if they do not terminate execution
warning('off', 'stats:kmeans:FailedToConverge');
%%  1. Signal Estimation
n_states = length(unique(idealization));
%   Idealization must find more than one state
%   in order to estimate signal
if n_states == 1
    disp("Idealization only found 1 state: SNR could not be estimated");
    return
end
%   Want to find all the state transitions in the
%   idealization
state_transitions = sort(abs(diff(idealization)));
state_transitions = state_transitions(state_transitions ~= 0);
%   Binary k-Means cluster of the state transition
%   amplitudes
[~, means] = kmeans(state_transitions, 2);
%   Signal estimate is the mean of the higher amplitude
%   cluster
signal = max(means);
%%  2. Noise Estimation
residuals = timeseries - idealization;
noise = std(residuals);
%%  3. SNR Calculation
estimated_SNR = signal / noise;
end
