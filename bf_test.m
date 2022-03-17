%% SIMO Beamforming Performance Test

rng(0); % Random, but repeatable!

% Generate detailed plots for this single simulation
[err_est, err_mea, old_err_meas] = calc_error_rate("BPSK", 10^( 2.0 /10), 4, true);

snrs = -5:0.1:3; % SNRs to simulate
numrxs = [1, 2, 3, 4]; % RX antenna numbers to simulate
colors = ["b", "r", "g", "m"]; % a list of colors for plots

%% Plot error rate vs SNR for different number of RX antennas
err_est = zeros(size(snrs));
err_mea = zeros(size(snrs));
old_err_meas = zeros(size(snrs));
plotnames = {};
figure;
for r=numrxs
    for ii=1:numel(snrs)
        [err_est(ii), err_mea(ii), old_err_meas(ii)] = calc_error_rate("BPSK", 10^(snrs(ii)/10), r, false);
    end
    plotnames{end+1} = sprintf("%d RX antennas with equalization (estimated)", r);
    plotnames{end+1} = sprintf("%d RX antennas with equalization", r);
    semilogy(snrs, err_est, colors(find(numrxs==r)), ...
             snrs, err_mea, colors(find(numrxs==r))+'x');
    hold on;
end
xlabel('E_0/N_0')
ylabel('Error Rate')
title('BPSK Symbol-Error-Rate vs. SNR')
grid on;
legend(plotnames, 'location', 'southwest');

%% Plot last set with and without equalization
figure;
semilogy(snrs, err_est, 'm', ...
         snrs, err_mea, 'mx',...
         snrs, old_err_meas, 'ro');
xlabel('E_0/N_0')
ylabel('Error Rate')
title('BPSK Symbol-Error-Rate vs. SNR')
grid on;
legend(sprintf("%d RX antennas with equalization (estimated)", r), ...
       sprintf("%d RX antennas with equalization", r), ...
       sprintf("%d RX antennas without equalization", r), ...
       'location', 'southwest');

function [estim, meas, orig_meas] = calc_error_rate(ctype, e0pern0, numrx, genPlot)
    if ctype=="BPSK"
        const = [-1, 1];        
    elseif ctype=="BFSK"
        const = [1i, 1];
    elseif ctype=="OOK"
        const = [0, sqrt(2)];
    else
        error("Const not supported");
    end

    d = abs(const(2) - const(1));
    estim = Q(d/sqrt(2) * sqrt(e0pern0*numrx));

    % Generate data
    N = ceil((3*10)^2/estim);
    data_bits = randi(2, 1,N);
    data = const(data_bits);

    % Generate Reference Signal
    N_RS = 400;
    rs_bits = repmat([1,2], 1, N_RS/2);
    rs = const(rs_bits);

    % Generate channel (phase shift per rx ant + noise + scaling)
    phase_shifts = exp(1i*(rand(numrx, 1)-0.5)*pi);
    rs_y   = repmat(  rs, numrx, 1) .* phase_shifts;
    data_y = repmat(data, numrx, 1) .* phase_shifts;    
    % Per channel scaling
    chan_scale = (rand(numrx, 1)*1.5 + 0.5);
    chan_scale = chan_scale/sqrt(mean(chan_scale.^2)); % Normalize
    %chan_scale = 0.5;
    rs_y   =   rs_y .* chan_scale;
    data_y = data_y .* chan_scale;
    n0 = 1/e0pern0;
    rs_y   =   rs_y + (sqrt(n0/2) * (randn(numrx,N_RS) + 1i*randn(numrx,N_RS)));
    data_y = data_y + (sqrt(n0/2) * (randn(numrx,N)    + 1i*randn(numrx,N)));    

    % Receiver processing
    x_eq = bf_equalizer(rs_y, data_y, rs);

    if ctype=="BPSK"
        err_idx = find( sign(real(x_eq)) ~= sign(real(data)) );
    elseif ctype=="BFSK"
        err_idx = find( sign(real(x_eq)-imag(x_eq)) ~= sign(real(data)-imag(data)) );
    elseif ctype=="OOK"
        err_idx = find( sign(real(x_eq)-sqrt(2)/2) ~= sign(real(data-sqrt(2)/2)) );
    else
        error("Const not supported");
    end
    errs = numel(err_idx);    
    meas = errs/N;
    
    if genPlot
        const_plot(data_y, x_eq, err_idx);
    end

    % Calculate errors without equalization
    if numrx >1
        x_mean = mean(data_y); % Average 'numrx' data points
    else
        x_mean = data_y;
    end
    
    if ctype=="BPSK"
        err_idx = find( sign(real(x_mean)) ~= sign(real(data)) );
    elseif ctype=="BFSK"
        err_idx = find( sign(real(x_x_meaneq)-imag(x_mean)) ~= sign(real(data)-imag(data)) );
    elseif ctype=="OOK"
        err_idx = find( sign(real(x_mean)-sqrt(2)/2) ~= sign(real(data-sqrt(2)/2)) );
    else
        error("Const not supported");
    end
    errs = numel(err_idx);    
    orig_meas = errs/N;
end

%% Plot constilation with noise
function const_plot(y, y_eq, err_idx)
   ysize = size(y);
   nrx = ysize(1);   
   t = tiledlayout(1,nrx);
   t.TileSpacing = 'compact';
   t.Padding = 'compact';
   for rx=1:nrx       
       nexttile;
       plot(real(y(rx,:)), imag(y(rx,:)), 'r.');
       grid on;
       ylim([-6, 6]); xlim([-6, 6]);
       title(sprintf('RX symbols antenna %d', rx))
   end   

   figure;
   plot(real(y_eq), imag(y_eq), '.');
   grid on; hold on;
   plot(real(y_eq(err_idx)), imag(y_eq(err_idx)), 'ro');    
   ylim([-3, 3]); xlim([-3, 3]);
   title('Equalized symbols with errors highlighted')
end

%% Gaussian Q function
function q = Q(x)
    q = 0.5*erfc(x/sqrt(2));
end