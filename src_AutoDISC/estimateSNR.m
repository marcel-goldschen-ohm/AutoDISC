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
%       1. Noise Estimation: Residuals between the
%       idealization and the noisy time series are
%       calculated, with their standard deviation
%       representing the average noise level of the trace
%       2. Signal Estimation: Step heights are dealt into an
%       array, weighted by the dwell time (in frames) preceding
%       and following them. All step heights within 2 noise
%       levels are discarded, and the mean of the reamining
%       step heights is taken as the signal estimate.
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

%%  1. Noise Estimation
residuals = timeseries - idealization;
noise = std(residuals);

%%  2. Signal Estimation
n_states = length(unique(idealization));
%   Idealization must find more than one state
%   in order to estimate signal
if n_states == 1
    disp("Idealization only found 1 state: SNR could not be estimated (Estimated SNR = 0)");
    estimated_SNR = 0;
    return
end
%   Want to find all the state transitions in the
%   idealization and the step amplitudes for each
state_transitions = abs(diff(idealization));
starts = unique([1; find(state_transitions ~= 0)]);
stops = [starts(2:end)-1; length(state_transitions)];
durs = (stops - starts) + 1;
step_heights = [];
for i = 1:numel(starts)
    if i == 1
        step_heights = [step_heights, repmat(state_transitions(starts(i)), 1, durs(i))];
    else
        step_heights = [step_heights, repmat(state_transitions(starts(i)), 1, durs(i)+durs(i-1))];
    end
end
step_heights = step_heights(step_heights > 2*noise);
signal = mean(step_heights);

%%  3. SNR Calculation
estimated_SNR = signal / noise;
if isnan(estimated_SNR)
    estimated_SNR = 0;
end
end
