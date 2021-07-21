function simulation = generateSMData(duration_f, n_sites, SNR, Q_f, state_amplitudes, state_het_dist, noise_profile)
%%  generateSMData Overview
%   Author: Argha Bandyopadhyay
%   Generates time series typical of single-molecule (SM) experiments
%   Inputs:
%       duration_f: number of points in final time series, 
%       in units of frames
%       n_sites: number of identical binding sites, each being 
%       simulated with given states/kinetics and then summed together to
%       form final output
%       SNR: signal-to-noise ratio (different definitions depending on
%       noise_profile, see below)
%       Q_f: transition rate matrix in units of frames^-1 - transitions
%       should be listed FROM row TO col and list 0 on diagonal
%       state_amplitudes: vector of the intensity values observed for each state at a
%       given site
%       state_het_dist: probability distribution describing 
%       the variation in intensity values for different
%       events in a given state (leave [] if no variation observed)
%       noise_profile: "gaussian" or "poisson"
%           if Gaussian, then SNR is defined such that the mean difference
%           in state amplitudes (step height) divided by the standard 
%           deviation of the Gaussian noise is equal to the SNR
%           if Poisson, then the SNR is defined such that the number of
%           photons observed at an arbitrary intensity level of 1 is equal
%           to (SNR^2)/0.95 (See Borner, et. al 2018 for derivation)
%   Output:
%       simulation: struct with following fields
%           ideal_data: time series of state amplitudes
%           noisy_data: time series of state amplitudes with applied noise
%           ideal_data_het: time series of state amplitudes with applied
%           per-event heterogeneity
%           noisy_data_het: time series of state amplitudes with applied
%           per-event heterogeneity and noise

%%  Implementation
% robustness/set up necessary parameters
duration_f = round(duration_f);
n_sites = round(n_sites);
n_states = size(Q_f, 1);
noise_profile = lower(noise_profile);

% unitary transition probability matrix
Q_f = Q_f - diag(sum(Q_f,2));
S = [Q_f ones(n_states,1)];
u = ones(1,n_states);
% equilibrium state probabilities
p0 = u*inv(S*(S'));


% dwell time dists and transition probability matrix
dwell_dists = cell(1, n_states);
trans_prob_mat = zeros(size(Q_f));
for i = 1:n_states
    rate_row = Q_f(i, :);
    % only want leaving rates (not from -> to same state)
    leaving_rates = rate_row(1:n_states ~= i);
    % dwell time distributions are defined as exponentials with
    % mean = 1 / sum(leaving rates)
    mu = 1 / sum(leaving_rates);
    dwell_dists{i} = makedist('Exponential', 'mu', mu);
    % transition probability matrix, with transition probability
    % to a state being the rate of leaving to the state / sum(all
    % leaving rates)
    trans_probs = leaving_rates / sum(leaving_rates);
    trans_prob_mat(i, 1:n_states ~= i) = trans_probs;
end

% pre-allocate space for outputs
ideal_data = zeros(1, duration_f);
ideal_data_het = zeros(1, duration_f);
ideal_data_mapped = zeros(1, duration_f);
noisy_data = zeros(1, duration_f);
noisy_data_het = zeros(1, duration_f);

% loop over number of sites
for site = 1:n_sites
    % simulate state sequence and corresponding dwell times
    % first, choose starting state and starting dwell
    draw = rand(1);
    event = 1;
    state_seq = [];
    dwell_seq = [];
    state_seq(event) = find(draw < cumsum(p0), 1);
    dwell_seq(event) = random(dwell_dists{state_seq(event)});
    % next, simulate as long as it is less than desired simulation
    % length
    while cumsum(dwell_seq) < duration_f
        cstate = state_seq(event);
        cstate_trans_probs = trans_prob_mat(cstate, :);
        trans_state_thresholds = cumsum(cstate_trans_probs);
        draw = rand(1);
        event = event + 1;
        state_seq(event) = find(draw < trans_state_thresholds, 1);
        dwell_seq(event) = random(dwell_dists{state_seq(event)});
    end
    % truncate last dwell if necessary
    cum_frames = cumsum(dwell_seq);
    extra_frames = cum_frames(end) - duration_f;
    dwell_seq(end) = dwell_seq(end) - extra_frames;
    
    % next, assign amplitudes for each event (with and without
    % state heterogeneity)
    n_events = length(dwell_seq);
    amplitude_seq = zeros(1, n_events);
    amplitude_seq_het = zeros(1, n_events);
    mean_step_amplitude = mean(diff(state_amplitudes));
    for st = 1:n_states
        state_mask = state_seq == st;
        n_st_events = sum(state_mask);
        amplitude_seq(state_mask) = state_amplitudes(st);
        het_changes = mean_step_amplitude*random(state_het_dist, [n_st_events,1])/100;
        het_direction = double(rand([n_st_events,1]) >= 0.5)*2 - 1;
        amplitude_seq_het(state_mask) = state_amplitudes(st) + het_changes.*het_direction;
    end
    
    % next, make idealized sequence based on amplitude + dwell
    % sequences
    ideal_seq = zeros(1, duration_f);
    ideal_seq_het = zeros(1, duration_f);
    ideal_seq_mapped = zeros(1,duration_f);
    frame_contribution = zeros(1, n_states);
    frame_counter = zeros(1, duration_f);
    frame_num = 1;
    % for each event, fill out the last frame, and then iterate
    % through the rest of the dwell time for the amplitude of interest
    for evt = 1:n_events
        frames_dur = dwell_seq(evt);
        while frames_dur > 0
            current_frame_remaining = 1 - frame_counter(frame_num);        
            % frame_frac will range from 0 to 1 to scale the
            % amplitude
            frame_frac = min(frames_dur, current_frame_remaining);
            % need to keep track of how much of the given frame has already
            % been used using frame_counter
            frame_counter(frame_num) = frame_counter(frame_num) + frame_frac;
            % ideal sequence values for the current frame are the amplitude
            % of the event scaled by frame_weight
            ideal_seq(frame_num) = ideal_seq(frame_num) + amplitude_seq(evt)*frame_frac;
            ideal_seq_het(frame_num) = ideal_seq_het(frame_num) + amplitude_seq_het(evt)*frame_frac;
            % keep track of which state is contributing to this event
            frame_contribution(state_seq(evt)) = frame_contribution(state_seq(evt)) + frame_frac;
            % the mapped sequence is the one with the most
            % contributing states for each frame
            [~, most_contributing_state] = max(frame_contribution);
            ideal_seq_mapped(frame_num) = state_amplitudes(most_contributing_state);
            % need to remove however many frames (from 0 to 1) we used
            frames_dur = frames_dur - frame_frac;
            % if the current frame has been filled completely, move to the
            % next frame
            if frame_counter(frame_num) == 1
                frame_contribution = zeros(1, n_states);
                frame_num = frame_num + 1;
                if frame_num > duration_f
                    break
                end
            end
        end
    end
    
    % next, apply noise based on selected noise profile to events that
    % aren't solely baseline (ideal ~= 0)
    noisy_seq = zeros(1, duration_f);
    noisy_seq_het = zeros(1, duration_f);
    no_BL_mask = ideal_seq ~= 0;
    if noise_profile == "gaussian"
        % Gaussian noise defined in relation to step heights
        signal = mean(diff(state_amplitudes));
        sigma_noise = signal / SNR;
        noisy_seq(no_BL_mask) = normrnd(ideal_seq(no_BL_mask), sigma_noise);
        noisy_seq_het(no_BL_mask) = normrnd(ideal_seq_het(no_BL_mask), sigma_noise);
    end
    
    % add sequence for current site to overall data
    ideal_data = ideal_data + ideal_seq;
    ideal_data_het = ideal_data_het + ideal_seq_het;
    ideal_data_mapped = ideal_data_mapped + ideal_seq_mapped;
    noisy_data = noisy_data + noisy_seq;
    noisy_data_het = noisy_data_het + noisy_seq_het;
end

% finally, apply noise to baseline (ideal == 0)
only_BL_mask = ideal_data == 0;
if noise_profile == "gaussian"
    % Gaussian noise defined in relation to step heights
    signal = mean(diff(state_amplitudes));
    sigma_noise = signal / SNR;
    noisy_data(only_BL_mask) = normrnd(ideal_data(only_BL_mask), sigma_noise);
    noisy_data_het(only_BL_mask) = normrnd(ideal_data_het(only_BL_mask), sigma_noise);
elseif noise_profile == "poisson"
    signal = mean(diff(state_amplitudes));
    baseline_photons = 50;
    baseline_noise = sqrt(baseline_photons);
    scaling_factor = SNR*baseline_noise/signal;
    ideal_data = baseline_photons + (ideal_data*scaling_factor);
    ideal_data_het = baseline_photons + (ideal_data_het*scaling_factor);
    ideal_data_mapped = baseline_photons + (ideal_data_mapped*scaling_factor);
    noisy_data = poissrnd(ideal_data);
    noisy_data_het = poissrnd(ideal_data_het);
end
simulation.ideal_data = ideal_data;
simulation.ideal_data_het = ideal_data_het;
simulation.ideal_data_mapped = ideal_data_mapped;
simulation.noisy_data = noisy_data;
simulation.noisy_data_het = noisy_data_het;
end