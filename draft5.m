close all;
clear all;
graphics_toolkit ("gnuplot");

# figure; hold on;
# t=0:0.01:1;
# ph=0;
# for i = 1:30
#   stem(i,mean(sin(2*pi*t).*sin(2*pi*t+ph)));
#   stem(i,mean(cos(2*pi*t).*sin(2*pi*t+ph)),'r');
#   ph+=0.1;
# endfor;
# return;

  %---------- Choice of parameters ------- %
  choice1 = menu ("Which number of bits in the bitstream do you want?", "12", "64", "102", "1020", "10200"); %number of bits in the FIFO
    if choice1 == 1
      n = 12;
    elseif choice1 == 2
      n = 64;
    elseif choice1 == 3
      n = 102;
    elseif choice1 == 4
      n = 1020;
    elseif choice1 == 5
      n = 10200;
    endif

  choice2 = menu ("Which QAM do you want to use?", "4", "16", "64");  %number of symbols in the constellation
    if choice2 == 1
      M = 4;
    elseif choice2 == 2
      M = 16;
    elseif choice2 == 3
      M = 64;
    endif

  bitsBySymbol = log2(M); %number of bits in one symbol

  %----------------- QAM Map Table --------------------

numberOfSymbols = n/bitsBySymbol; %number of symbols in the FIFO

% Define mapping table applying Gray mapping

if M == 4
  mappingTable(1) =  1 + 1*j;
  mappingTable(2) = -1 + 1*j;
  mappingTable(3) =  1 - 1*j;
  mappingTable(4) = -1 - 1*j;

elseif M == 16
  mappingTable(1:4) = -3;
  mappingTable(5:8) = -1;
  mappingTable(9:12) = +3;
  mappingTable(13:16) = +1;
  for i = 0:15
    if mod(i,4) == 0
      mappingTable(i+1) = mappingTable(i+1) -3*j;
    elseif mod(i+3,4) == 0
      mappingTable(i+1) = mappingTable(i+1) -1*j;
    elseif mod(i+1,4) == 0
      mappingTable(i+1) = mappingTable(i+1) +1*j;
    elseif mod(i+2,8) == 0
      mappingTable(i+1) = mappingTable(i+1) +3*j;
    endif
  endfor


elseif M == 64
  mappingTable(1:8)   = + 7*j;
  mappingTable(9:16)  = + 5*j;
  mappingTable(17:24) = + 1*j;
  mappingTable(25:32) = + 3*j;
  mappingTable(33:40) = - 7*j;
  mappingTable(41:48) = - 5*j;
  mappingTable(49:56) = - 1*j;
  mappingTable(57:64) = - 3*j;

  for i = 0:63
    if mod(i+2,8) == 0
      mappingTable(i+1) = mappingTable(i+1) +1;
    elseif mod(i+1,8) == 0
      mappingTable(i+1) = mappingTable(i+1) +3;
    elseif mod(i+3,8) == 0
      mappingTable(i+1) = mappingTable(i+1) +5;
    elseif mod(i+4,8) == 0
      mappingTable(i+1) = mappingTable(i+1) +7;
    elseif mod(i+6,8) == 0
      mappingTable(i+1) = mappingTable(i+1) -1;
    elseif mod(i+5,8) == 0
      mappingTable(i+1) = mappingTable(i+1) -3;
    elseif mod(i+7,8) == 0
      mappingTable(i+1) = mappingTable(i+1) -5;
    elseif mod(i,8) == 0
      mappingTable(i+1) = mappingTable(i+1) -7;
    endif
  endfor
endif



% ---- ---- Simulation Paramters ---- ---- %
k=1;
end_k = numberOfSymbols +1;

dt = 1e-8;
end_t = 1e-3;

window_size = 5;

PLOT_TX = 0;
PLOT_RX = 0;

dataIn = randi(2,n,1) - 1; % Generate vector of binary data
dataOut = zeros(n,1);

if PLOT_TX
  figure(1);
  hold on;
endif;
if PLOT_RX
  figure(2);
  hold on;
endif;
% ---- ---- Simulation Paramters ---- ---- %


% ---- ---- Transceiver Paramters ---- ---- %
fc               = 1e6;
# fc_tx            = fc*1.0001;
# fc_tx            = fc*1.00002; % typical error 20ppm
fc_tx            = fc;
ph_tx            = 0.4135;
fc_rx            = fc;
ph_rx            = 0.0;
# Tsymb            = 1e-4;     % Symbol period
Tsymb            = 10/fc;     % Symbol period = 100/fc = 100 fc periods
samples_per_symb = 16;
Ts               = Tsymb/samples_per_symb; % Sampling period: 16 samples/symbol
alpha_param      = 1;
beta_param       = 0.25;
delay            = Tsymb/3;
delay            = 0.0;
SNR              = 100;
  % ---- ---- Pulse shaping ---- ---- %
  t = 0:dt:Tsymb-dt;
  tl= length(t);
#   pulse_shaping = sin(pi*(0:dt:Tsymb-dt)/(Tsymb-dt));
  pulse_shaping = 0.5+sin(pi*(0:dt:Tsymb-dt)/(Tsymb-dt)) + 0.5;
#   pulse_shaping = ones(1,tl);
# SNR = 100;
% ---- ---- Transceiver Paramters ---- ---- %






# % ---- ---- Transceiver Paramters ---- ---- %
#   % ---- ---- Bad Paramters ---- ---- %
#   delay = 0.005e-3;
#   ph_tx = 0;
#   ph_rx = 2*pi*0.13;
#   ph_tx = 0.451;
#   fc_tx = fc;
#   fc_rx = fc*1.01;
#   SNR = 0;
#   % ---- ---- Good Paramters ---- ---- %
# #   delay = 0.0;
# #   ph_tx = 0;
#   ph_rx = 0;
#   fc_tx = fc;
# #   fc_rx = fc;
#   SNR = 100;
#
#
#
#   maf_f_l = 100;
#   fifo_f = zeros(1,maf_f_l);
#   maf_f(1) = 0;
#
#   % ---- ---- Pulse shaping ---- ---- %
#   pulse_shaping = sin(2*pi*(0:dt:Tsymb-dt)/(Tsymb-dt));
#   pulse_shaping = ones(1,tl);
#
% ---- ---- END Transceiver Paramters ---- ---- %





xi_aux_t = zeros(1,tl);
xq_aux_t = zeros(1,tl);

x2_t = zeros(1,2*tl);
y3_t = zeros(1,3*tl);
yi2_t = zeros(1,2*tl);
yq2_t = zeros(1,2*tl);

aux_error = 0;
last_phase = 0;
sync_flag = -20;
freeze = 1;

sync_sample = 0;

% ---- ---- Main Loop ---- ---- %
while(k < end_k)

  % ---- ---- Transmitter ---- ---- %

  if sync_flag <= 0 %Dummy symbol to lock the PLL
    xi_k = 1;
    xq_k = 1;
  else
    symbolBits = dataIn((k-1)*bitsBySymbol+1:(k-1)*bitsBySymbol+bitsBySymbol);

    if M == 4
      symbolIndex = 2^1 * symbolBits(1) + 2^0 * symbolBits(2);
    elseif M == 16
      symbolIndex = 2^3 * symbolBits(1) + 2^2 * symbolBits(2) + 2^1 * symbolBits(3) + 2^0 * symbolBits(4);
    elseif M == 64
      symbolIndex = 2^5 * symbolBits(1) + 2^4 * symbolBits(2) + 2^3 * symbolBits(3) + 2^2 * symbolBits(4) +  2^1 * symbolBits(5) + 2^0 * symbolBits(6);
    endif
    symbolIndex;
     % Mapping
    symbol = mappingTable(symbolIndex + 1);
    xi_k = real(symbol);
    xq_k = imag(symbol);
  endif

  xsi_t = pulse_shaping .* xi_k;
  xsq_t = pulse_shaping .* xq_k;

  xi_t = xsi_t .* cos(2*pi*fc_tx*t+ph_tx);
  xq_t = xsq_t .* sin(2*pi*fc_tx*t+ph_tx);

  x_t = xi_t + xq_t;
  x2_t = [x2_t(1+tl:2*tl) x_t];
  % ---- ---- END Transmitter ---- ---- %


  % ---- ---- Receiver ---- ---- %;

  %CHANNEL
  delay_idx = mod(tl-floor(delay/dt),tl)+1;
  idx = [delay_idx:(delay_idx+tl-1)];
  y3_t = [y3_t(1+tl:3*tl) awgn(x2_t(idx), SNR)];
  y2_t = y3_t(1:2*tl);
#   y_t = awgn(x_t, SNR);
#   y_t = awgn(x2_t(idx), SNR);
#   y_t = awgn(x2_t, SNR);

  %CARRIER FREQUENCY AND PHASE RECOVERY

#   carrier_inph = cos(2*pi*fc_rx*t + ph_rx);
#   carrier_quad = sin(2*pi*fc_rx*t + ph_rx);

#   yi_t = y_t .* carrier_inph;
#   yq_t = y_t .* carrier_quad;
#   % FILTERING BY FFT
# #   aux = yi_t + j*yq_t;
# #   aux = fft(aux) .* [1 1 zeros(1,tl-2)];
# #   aux = ifft(aux);
# #   y_t = aux;
#   yi2_t = [yi2_t(1+tl:2*tl) yi_t];
#   yq2_t = [yq2_t(1+tl:2*tl) yq_t];
  t2 = t(1):dt:(t(end)+Tsymb);
  carrier_inph = cos(2*pi*fc_rx*t2 + ph_rx);
  carrier_quad = sin(2*pi*fc_rx*t2 + ph_rx);
  yi2_t = y2_t .* carrier_inph;
  yq2_t = y2_t .* carrier_quad;


#   sample_number = mod(floor(1/fc/4/dt) + floor(ph_rx/2/pi*Tsymb/dt), Tsymb/dt);
#   sample_number = mod(floor(1/fc/4/dt), Tsymb/dt);
  indx = Ts*sync_sample:Ts:(Tsymb-Ts)+Ts*sync_sample;
  indx /= dt;
  indx = floor(indx);
  indx += 1;
  yi_s = yi2_t(indx);
  yq_s = yq2_t(indx);

  indx = floor(mod(indx-1, Tsymb/dt)+1) ;

  yi_k = mean(yi_s./pulse_shaping(indx));
  yq_k = mean(yq_s./pulse_shaping(indx));

  if sync_flag <= 0
    sync_flag ++;
    freeze = 1;
    symbolIndexAfter = 0;
  else
    receivedSymbols = yi_k + j*yq_k;
    [mindiff minIndex] = min(receivedSymbols - mappingTable);
    symbolIndexAfter = minIndex - 1;

    aux_symbolIndexAfter = symbolIndexAfter;
    for i = 1:bitsBySymbol
      bits(i) = mod(aux_symbolIndexAfter,2);
      aux_symbolIndexAfter = floor(aux_symbolIndexAfter/2);
    end
    bits = fliplr(bits);
    for i = 1:bitsBySymbol
      dataOut((k-1)*bitsBySymbol +i) = bits(i);
    end;
    sync_sample = mod(sync_sample, samples_per_symb);
  endif;

  mappingTable(symbolIndexAfter+1)
  yi_k + j * yq_k
  arg(mappingTable(symbolIndexAfter+1))
  arg(yi_k + j * yq_k)
  phase_error = arg(mappingTable(symbolIndexAfter+1)) - arg(yi_k + j * yq_k)
  aux_error += phase_error;
  phase_correction = alpha_param * phase_error + beta_param * aux_error;

  ph_rx = mod(ph_rx + phase_correction, 2*pi)


    figure(10); hold on;
#     plot(t,aux_phase_error_t);
    plot([t(1) t(1)+Tsymb],phase_error*[1 1],'k');
    plot([t(1) t(1)+Tsymb],ph_rx*[1 1],'g');
#     plot([t(1) t(1)+Tsymb],fc_error*[1 1],'r');
%     plot(t,fc_rx*ones(1,length(t)),'c');
    plot([t(1) t(1)+Tsymb],aux_error*[1 1],'m','linewidth',2);
    plot([t(1) t(1)+Tsymb],phase_correction*[1 1],'c','linewidth',2);


    % ---- ---- END Receiver ---- ---- %




    % ---- ---- Plot Transmitter Signals ---- ---- %
    if (PLOT_TX == 1)
      figure(1);
#       subplot(4,1,1); hold on;
#       plot([t(1) t(1)+Tsymb],[0,0],'k-');
      stem(t(1), xi_k,'b','linewidth',3);
      stem(t(1), xq_k,'r','linewidth',2);
      axis([t(1)-window_size*Tsymb t(1)+Tsymb -1.5 1.5]); %TODO: implement a plotting buffer with a windown from 5 to 10 Tsymb

#       subplot(4,1,2);  hold on;
      plot(t,xsi_t-3,'b:','linewidth',2);
      plot(t,xsq_t-3,'r-','linewidth',1);
      axis([t(1)-window_size*Tsymb t(1)+Tsymb -4.5 1.5]); %TODO: implement a plotting buffer with a windown from 5 to 10 Tsymb

#       subplot(4,1,3);  hold on;
#       plot(t,xi_t-6,'b:','linewidth',1);
#       plot(t,xq_t-6,'r-','linewidth',1);
#       axis([t(1)-window_size*Tsymb t(1)+Tsymb -7.5 1.5]); %TODO: implement a plotting buffer with a windown from 5 to 10 Tsymb

#       subplot(4,1,4); hold on;
#       plot(t,x_t,'g-','linewidth',1);
#       axis([t(1)-window_size*Tsymb t(1)+Tsymb -1.5 1.5]); %TODO: implement a plotting buffer with a windown from 5 to 10 Tsymb

        drawnow ("expose");
      title('xik(blue)/xqk(red), xsit/q, xit/q, xt');
    endif;
    % ---- ---- END Plot Transmitter Signals ---- ---- %


    % ---- ---- Plot Receiver Signals ---- ---- %
    if (PLOT_RX == 1)
      figure(2); hold on;
#       subplot(4,1,1); hold on;
#       plot(t,y_t,'g-','linewidth',1)
# #       plot([t(1) t(1)+Tsymb],[0,0],'k-');

#       subplot(4,1,2); hold on;
#       plot(t,yi_t-6,'b:','linewidth',1)
#       plot(t,yq_t-6,'r-','linewidth',1)

#       subplot(4,1,3); hold on;
#       stem(t(indx),yi_s,'m:','linewidth',2)
#       stem(t(indx),yq_s,'g-','linewidth',1)
      plot(t(indx),yi_s,'m:','linewidth',2)
      plot(t(indx),yq_s,'g-','linewidth',1)
      axis([t(1)-window_size*Tsymb t(1)+Tsymb -1.5 1.5]); %TODO: implement a plotting buffer with a windown from 5 to 10 Tsymb

#       figure(3); hold on;
#       subplot(4,1,4); hold on;
      stem(t(1), yi_k,'bx','linewidth',3,'linewidth',4,'markersize',15);
      stem(t(1), yq_k,'rx','linewidth',2,'linewidth',4,'markersize',15);
#       stem(t(indx), yi_s./pulse_shaping(indx),'bx','linewidth',3,'linewidth',4,'markersize',15);
#       stem(t(indx), yq_s./pulse_shaping(indx),'rx','linewidth',2,'linewidth',4,'markersize',15);
      axis([t(1)-window_size*Tsymb t(1)+Tsymb -1.5 1.5]); %TODO: implement a plotting buffer with a windown from 5 to 10 Tsymb

#       axis([t(1)-window_size*Tsymb t(1)+Tsymb -11 2]); %TODO: implement a plotting buffer with a windown from 5 to 10 Tsymb
#       title('yik(blue)/yqk(red), yifiltert/q, yit/q, yt');
      drawnow ("expose");
    endif;
    % ---- ---- END Plot Receiver Signals ---- ---- %

    % ---- ---- Update time and plots ---- ---- %
    t = t + Tsymb;
    if freeze == 0
      k++;
    elseif sync_flag == 1
      freeze = 0;
    end
    fflush(stdout);
    % ---- ---- END Update time and plots ---- ---- %
endwhile;
% ---- ---- END Main Loop ---- ---- %

dataInTest = dataIn(1:end-6);
dataOutTest = dataOut(7:end);

%-------------- BER -------------

total_error = 0;

for i = 1:n-6
    if dataIn(i) != dataOut(i+6)
      total_error ++;
    end
end

% Calculation of BER to return the result
BER = total_error/n;

% Showing final results
disp(['Total wrong bits = ' num2str(total_error)]);
disp(['BER = ' num2str(BER)]);
%
%figure(3);
%plot(vco_f,'c;vco f;');
