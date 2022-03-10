

[err_est, err_mea] = calc_error_rate("BPSK", 10^( 0.0 /10), 4);
fprintf('orig err est = %d\n', err_est);
fprintf('measured error = %d\n', err_mea);


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

    % Generate data
    N = ceil((3*10)^2/estim);
    data_bits = randi(2, 1,N);
    data = const(data_bits);

    % Generate Reference Signal
    N_RS = 400; % TODO try different nums
    rs_bits = repmat([1,2], 1, N_RS/2);
    rs = const(rs_bits);

    % Generate channel (phase shift per rx ant + noise + scaling)
    phase_shifts = exp(1i*(rand(numrx, 1)-0.5)*pi);
    rs_y   = repmat(  rs, numrx, 1) .* phase_shifts;
    data_y = repmat(data, numrx, 1) .* phase_shifts;
    n0 = 1/e0pern0;
    rs_y   =   rs_y + (sqrt(n0/2) * (randn(numrx,N_RS) + 1i*randn(numrx,N_RS)));
    data_y = data_y + (sqrt(n0/2) * (randn(numrx,N)    + 1i*randn(numrx,N)));
    % Per channel scaling
    chan_scale = (rand(numrx, 1)*1.5 + 0.5);
    rs_y   =   rs_y .* chan_scale;
    data_y = data_y .* chan_scale;
    % TODO add per channel noise power?

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
    
    const_plot(data_y, x_eq, err_idx);
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