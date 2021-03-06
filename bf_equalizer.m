function [x_data] = bf_equalizer(y_rs, y_data, x_rs)
    nrx = numel(y_rs)/numel(x_rs); % Num RX antennas
    %% Generate channel estimates
    h = y_rs .* conj(x_rs); % Channel estimate for each symbol
    H = mean(h,2); % Average to create a total channel estimate for each rx antenna    

    %% Generate noise estimates
    n = y_rs - H*x_rs; % Noise estimate for each symbol
    n_est = mean(abs(n),2); % Averate to get total noise estimate for each rx antenna

    %% Generate weights
    P = 1; % Assume signal power
    n_pwr = n_est'*n_est; % Noise power
    w = H' / (H*H' + n_pwr/P*eye(nrx)); % Minimized-Mean-Square-Error equasion
    % w = H' / (H*H' + eye(nrx)); % MMSE equation without SNR compensation

    %% Apply weights
    x_data = w*y_data; % Weights "undo" the channel H with minimal error
end

