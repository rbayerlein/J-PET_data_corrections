function [t_start,t_len,t_scaled,t_acq] = GetFramingInfo(t_vec,t12)
%GETFRAMINGINFO Summary of this function goes here
%   Detailed explanation goes here

if nargin ~= 2 || mod(numel(t_vec),2) ~= 0
    error('Usage: [t_start,t_len,t_acq] = GetFramingInfo([n_frames,frame_length,[n_frames,frame_length,...]])');
end

t_start = []; % vector of frame start times (in s) - initialization
t_len = []; % vector of frame lengths (in s) - initialization
t_scaled = []; % vector of t12-scaled simulation time (in s) - initialization

t_now = 0; % current time (in s) - initialization
for i = 1:2:numel(t_vec) % get n_frames
    n_frames = t_vec(i);
    for j = 1:n_frames
        t_start(end+1) = t_now; % C++11 style vector::push_back
        t_len(end+1) = t_vec(i+1); % C++11 style vector::push_back
        % scale with effective decay factor
        % see Cherry (pg. 36) Eqn 4.15
        x = log(2)*t_vec(i+1)/t12; % Eqn 4.16
        % Based on eqn 4.15: ad/a0 = delta_t_prime/delta_t
        % => delta_t_prime = delta_t*ad/a0
        % => delta_t_prime = delta_t*DF_eff
        t_scaled(end+1) = t_vec(i+1)*(1-exp(-x))/x; % C++11 style vector::push_back
        t_now = t_now + t_vec(i+1);
    end
end

t_acq = sum(t_vec(1:2:end).*t_vec(2:2:end)); % total acquisition time (in s)

end

