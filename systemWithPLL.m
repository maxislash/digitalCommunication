close all;
clear all;
graphics_toolkit ("gnuplot");

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
last_phase = 0;
sync_flag = -60;
freeze = 1;
i = 2;
k=1;
end_k = numberOfSymbols +1;
Ts = 1e-3;
dt = 1e-5;
end_t = 1e-3;
t = 0:dt:Ts-dt;
lt = length(t);
window_size = 5;
PLOT_TX = 1 ;
PLOT_RX = 1;

% ---- ---- Transceiver Paramters ---- ---- %
fc = 1e3;
fs = 1/dt;
dataIn = randi(2,n,1) - 1; % Generate vector of binary data
dataOut = zeros(n,1);
  % ---- ---- Bad Paramters ---- ---- %
  delay = 0.005e-3;
  phase_tx = 0;
  phase_rx = 2*pi*0.13;
  phase_tx = 0.451;
  fc_tx = fc;
  fc_rx = fc*1.01;
  SNR = 0;
  % ---- ---- Good Paramters ---- ---- %
  delay = 0.0;
 % phase_rx = 0;
  %phase_tx = 0;
  fc_tx = fc;
 % fc_rx = fc;
  SNR = 100;


  maf_f_l = 100;
  fifo_f = zeros(1,maf_f_l);
  maf_f(1) = 0;

  % ---- ---- Pulse shaping ---- ---- %
  pulse_shaping = sin(2*pi*(0:dt:Ts-dt)/(Ts-dt));
  pulse_shaping = ones(1,lt);

% ---- ---- END Transceiver Paramters ---- ---- %



%if PLOT_TX
%  figure(1);
%  hold on;
%endif;
%if PLOT_RX
%  figure(2);
%  hold on;
%endif;


# xi_aux_t = zeros(1,floor(Ts/dt));
# xq_aux_t = zeros(1,floor(Ts/dt));
xi_aux_t = zeros(1,lt);
xq_aux_t = zeros(1,lt);

x2_t = zeros(1,2*lt);

% ---- ---- Main Loop ---- ---- %
while(k < end_k)


  % ---- ---- Transmitter ---- ---- %
#   xi_k = -xi_k;
#   xq_k = 1;

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

  xi_t = xsi_t .* cos(2*pi*fc_tx*t+phase_tx);
  xq_t = xsq_t .* sin(2*pi*fc_tx*t+phase_tx);

  x_t = xi_t + xq_t;
  x2_t = [x2_t(1+lt:2*lt) x_t];
  % ---- ---- END Transmitter ---- ---- %


  % ---- ---- Receiver ---- ---- %;

  %CHANNEL
  y_t = awgn(x_t, SNR);
  delay_idx = mod(floor(delay/dt),lt)+1;
  idx = [delay_idx:(delay_idx+lt-1)];
  y_t = awgn(x2_t(idx), SNR);

  %CARRIER FREQUENCY AND PHASE RECOVERY
#   if k == 0  %TODO : send the carrier only at the beginning to recover the frequency and the phase

    k;
    aux_vco_t = -sin(2*pi*fc_rx*t + phase_rx - pi/4);
    aux_phase_error_t = y_t .* aux_vco_t;
    phase_error = mean(aux_phase_error_t);

    fc_error = (phase_rx - last_phase);
    last_phase = phase_rx;
    
    maf_f(i) = maf_f(i-1) + fc_error/Ts; %Derivate the phase to have the frequency
    maf_f(i) -= fifo_f(maf_f_l);
    fifo_f(2:maf_f_l) = fifo_f(1:maf_f_l-1);
    fifo_f(1) = fc_error/Ts;

    if sync_flag  <= 0  %TODO : send the carrier only at the beginning to recover the frequency and the phase
      freeze = 1;
      phase_rx = phase_rx + phase_error;
      fc_rx = fc_rx + maf_f(i)/pi/2/maf_f_l;
  %     fc_rx = fc_rx + fc_error * 10   
      if phase_error < 0.001
        sync_flag ++;     
      endif;
      
    endif;

    i++;
    fc_rx;
    phase_rx;
%    figure(10); hold on;
%   % plot(t,aux_phase_error_t);
%    plot([t(1) t(1)+Ts],phase_error*[1 1],'k');
%    plot([t(1) t(1)+Ts],phase_rx*[1 1],'g');
%    plot([t(1) t(1)+Ts],fc_error*[1 1],'r');
%#     plot(t,fc_rx*ones(1,length(t)),'c');




#
#
#
#     delta_t = 1e-5;
#     max_t = 5e-1;
#     t_sync = 0:delta_t:max_t;
#     y_t = cos(2*pi*fc_rx*t_sync + phase_rx);
#
#     for i = 2:length(t_sync)
#
#       vco_t(i) = -sin(2*pi*vco_f(i-1)*t_sync(i) + vco_p(i-1));
#       aux_t(i) = y_t(i) * vco_t(i);
#       maf(i) = maf(i-1) + aux_t(i);
#       maf(i) -= fifo(maf_l);
#
#       fifo(2:maf_l) = fifo(1:maf_l-1);
#       fifo(1) = aux_t(i);
#
#       vco_p(i) = vco_p(i-1) + maf(i-1)/maf_l/100;
#
#       maf_f(i) = maf_f(i-1) + (vco_p(i) - vco_p(i-1))/delta_t; %Derivate the phase to have the frequency
#       maf_f(i) -= fifo_f(maf_f_l);
#       fifo_f(2:maf_f_l) = fifo_f(1:maf_f_l-1);
#       fifo_f(1) = (vco_p(i) - vco_p(i-1))/delta_t;
#
#       if (mod(i,10000) == 0)
#         vco_f(i) = vco_f(i-1) + maf_f(i)/pi/2/maf_f_l;
#       else
#         vco_f(i) = vco_f(i-1);
#       endif;
#
#   %    carrier_inph(step-1) = cos(2*pi*vco_f(step-1)*step*dt+vco_p(step-1));
#   %    carrier_quad(step-1) = sin(2*pi*vco_f(step-1)*step*dt+vco_p(step-1));
#     endfor
#     k++;
#
#     p_step = 100;
#
#     figure(4); hold on;
#     plot(t_sync(1:p_step:end),aux_t(1:p_step:end));
#     plot(t_sync(1:p_step:end),maf(1:p_step:end)/maf_l*2,'r;phase error;');
#     plot(t_sync(1:p_step:end),vco_p(1:p_step:end),'g;vco phase;');
#      %plot(t(1:p_step:end),vco_f(1:p_step:end)-rx_f(1:p_step:end),'c;vco f - rx f;')
#     plot(t_sync(1:p_step:end),maf_f(1:p_step:end)/1000/pi,'k;freq error;')
#
#     figure(5); hold on;
#     plot(t_sync(1:p_step:end),fc_rx*ones(1,length(t_sync(1:p_step:end))),'b;rx f;');
#     plot(t_sync(1:p_step:end),vco_f(1:p_step:end),'c;vco f;');

#   else

#    carrier_inph = cos(2*pi*vco_f(end)*t+vco_p(end));
#    carrier_quad = sin(2*pi*vco_f(end)*t+vco_p(end));

    %RECOVERY OF YI_T AND YQ_T
    carrier_inph = cos(2*pi*fc_rx*t + phase_rx);
    carrier_quad = sin(2*pi*fc_rx*t + phase_rx);

    yi_t = y_t .* carrier_inph;
    yq_t = y_t .* carrier_quad;
    % TODO The LBP filter is missing
  %   lbp_length = 4;
  %   lbp_filter = [ones(1,lbp_length) zeros(1,lt-lbp_length)];
  %   lbp_filter /= sum(lbp_filter);
  %   yi_t = filter(lbp_filter, [1 zeros(1,19)], [xi_aux_t yi_non_lbp_t])(lt+1:2*lt);
  %   yq_t = filter(lbp_filter, [1 zeros(1,19)], [xi_aux_t yq_non_lbp_t])(lt+1:2*lt);
  %   xi_aux_t = yi_non_lbp_t;
  %   xq_aux_t = yq_non_lbp_t;

  % FILTERING BY FFT
    aux = yi_t + j*yq_t;
    aux = fft(aux) .* [1 1 zeros(1,lt-2)];
    aux = ifft(aux);
    yi_filter_t = 2*real(aux);
    yq_filter_t = 2*imag(aux);
    yi_k = yi_filter_t(1);
    yq_k = yq_filter_t(1);

  % TODO : FIR filter but for the mean time, FFT is good to implement the carrier recovery
  %  %figure(3)
  %  b = fir1(100,0.001);
  %  %freqz(b,10e6);
  %  yi_filter_t = 2*filter(b,1,yi_t);
  %  yq_filter_t = 2*filter(b,1,yq_t);
  %  yi_k = yi_filter_t(50)
  %  yq_k = yq_filter_t(50);

  %QAM DEMAPPER
    receivedSymbols = yi_k + j*yq_k;
    [mindiff minIndex] = min(receivedSymbols - mappingTable);
    symbolIndexAfter = minIndex - 1;
   
    if sync_flag == 1
      if symbolIndexAfter == 0 
        'Phase test'
        if phase_error >= 0.001
          'Correction of the phase'
          phase_rx = mod(phase_rx + phase_error,2*pi);    
        endif
      endif 

      for i = 1:bitsBySymbol
        bits(i) = mod(symbolIndexAfter,2);
        symbolIndexAfter = floor(symbolIndexAfter/2);
      endfor
      bits = fliplr(bits);

      for i = 1:bitsBySymbol
        dataOut((k-1)*bitsBySymbol +i) = bits(i);
      endfor
    
    endif;
  


    % ---- ---- END Receiver ---- ---- %




%    % ---- ---- Plot Transmitter Signals ---- ---- %
%    if (PLOT_TX == 1)
%      figure(1);
%      plot([t(1) t(1)+Ts],[0,0],'k-');
%      stem(t(1), xi_k,'b','linewidth',3);
%      stem(t(1), xq_k,'r','linewidth',2);
%
%      plot([t(1) t(1)+Ts],[-3,-3],'k-');
%      plot(t,xsi_t-3,'b:','linewidth',2);
%      plot(t,xsq_t-3,'r-','linewidth',1);
%
%      plot([t(1) t(1)+Ts],[-6,-6],'k-');
%      plot(t,xi_t-6,'b:','linewidth',1);
%      plot(t,xq_t-6,'r-','linewidth',1);
%
%      plot([t(1) t(1)+Ts],[-9,-9],'k-');
%      plot(t,x_t-9,'g-','linewidth',1);
%      figure(1);
%      axis([t(1)-window_size*Ts t(1)+Ts -11 2]); %TODO: implement a plotting buffer with a windown from 5 to 10 Ts
%      drawnow ("expose");
%      title('xik(blue)/xqk(red), xsit/q, xit/q, xt');
%    endif;
%    % ---- ---- END Plot Transmitter Signals ---- ---- %
%
%
%    % ---- ---- Plot Receiver Signals ---- ---- %
%    if (PLOT_RX == 1)
%      figure(2);
%      plot(t,y_t-9,'g-','linewidth',1)
%      plot([t(1) t(1)+Ts],[-9,-9],'k-');
%      axis([t(1)-window_size*Ts t(1)+Ts -11 2]); %TODO: implement a plotting buffer with a windown from 5 to 10 Ts
%      drawnow ("expose");
%
%      plot([t(1) t(1)+Ts],[-6,-6],'k-');
%      plot(t,yi_t-6,'b:','linewidth',1)
%      plot(t,yq_t-6,'r-','linewidth',1)
%
%      plot([t(1) t(1)+Ts],[-3,-3],'k-');
%      plot(t,yi_filter_t-3,'b:','linewidth',2)
%      plot(t,yq_filter_t-3,'r-','linewidth',1)
%
%      plot([t(1) t(1)+Ts],[0,0],'k-');
%      stem(t(1), yi_k,'b','linewidth',3);
%      stem(t(1), yq_k,'r','linewidth',2);
%
%      title('yik(blue)/yqk(red), yifiltert/q, yit/q, yt');
%    endif;
    % ---- ---- END Plot Receiver Signals ---- ---- %

    % ---- ---- Update time and plots ---- ---- %
    t = t + Ts;
    if freeze == 0
      k++;
    elseif sync_flag == 1
      freeze = 0; 
    endif;
    fflush(stdout);
    % ---- ---- END Update time and plots ---- ---- %
#   endif
#   vco_p(end)
#   k
endwhile;
% ---- ---- END Main Loop ---- ---- %
dataInTest = dataIn(1:end-2);
dataOutTest = dataOut(3:end);

%-------------- BER -------------

total_error = 0;

for i = 1:n-2
    if dataIn(i) != dataOut(i+2)
      total_error ++;
    endif
endfor;

% Calculation of BER to return the result
BER = total_error/n;

% Showing final results
disp(['Total wrong bits = ' num2str(total_error)]);
disp(['BER = ' num2str(BER)]);
%
%figure(3);
%plot(vco_f,'c;vco f;');
