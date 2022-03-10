%close all;

%% Run simulations for the 3 constilations at different SNRs
snrs = -5:0.1:3;
numrxs = [1, 2, 3, 4];
colors = ["b", "r", "g", "m"];

err_est = zeros(size(snrs));
err_mea = zeros(size(snrs));
plotnames = {};
figure;
% hold on;
for r=numrxs
    for ii=1:numel(snrs)
        [err_est(ii), err_mea(ii)] = calc_error_rate("BPSK", 10^(snrs(ii)/10), r);
    end
    plotnames{end+1} = "";
    plotnames{end+1} = sprintf("%d RX antennas", r);    
    semilogy(snrs, err_est, colors(find(numrxs==r)), ...
             snrs, err_mea, colors(find(numrxs==r))+'x');
    hold on;
end
xlabel('E_0/N_0')
ylabel('Error Rate')
title('BPSK Symbol-Error-Rate vs. SNR')
grid on;
legend(plotnames, 'location', 'southwest');

%% Estimate the error rate, and perform an MC simulation
function [estim, meas] = calc_error_rate(ctype, e0pern0, numrx)
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

    N = ceil((3*10)^2/estim);
    data_bits = randi(2, 1,N);
    data = const(data_bits);
    
    % Generate and add AWGN
    n0 = 1/e0pern0;
    n = sqrt(n0/2) * (randn(numrx,N) + 1i*randn(numrx,N));
    y = repmat(data, numrx, 1) + n;
    
    % Receiver processing
    if numrx >1
        z = mean(y); % Average 'numrx' data points
    else
        z = y;
    end
    
    if ctype=="BPSK"
        err_idx = find( sign(real(z)) ~= sign(real(data)) );
    elseif ctype=="BFSK"
        err_idx = find( sign(real(z)-imag(z)) ~= sign(real(data)-imag(data)) );
    elseif ctype=="OOK"
        err_idx = find( sign(real(z)-sqrt(2)/2) ~= sign(real(data-sqrt(2)/2)) );
    else
        error("Const not supported");
    end
    errs = numel(err_idx);    
    meas = errs/N;

    % const_plot(y, err_idx);
end


%% Gaussian Q function
function q = Q(x)
    q = 0.5*erfc(x/sqrt(2));
end

%% Plot constilation with noise
function const_plot(y, err_idx)
   close all; figure;
   plot(real(y), imag(y), '.');
   grid on; hold on;
   plot(real(y(err_idx)), imag(y(err_idx)), 'ro'); 
   ylim([-2, 2]); xlim([-2, 2])
end