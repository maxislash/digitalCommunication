close all;
clear all;
graphics_toolkit ("gnuplot");

% ---- ---- Simulation Paramters ---- ---- %
k=0;
end_k = 8;
Ts = 1e-3;
dt = 1e-6;
end_t = 1e-3;
t = 0:dt:Ts-dt;
lt = length(t);
window_size = 8;
PLOT_TX = 1;
PLOT_RX = 1 ;


% ---- ---- Transceiver Paramters ---- ---- %
fc = 10e3;
i_seq = [1 -1 -1 1];
q_seq = [1 1 -1 -1];
  % ---- ---- Bad Paramters ---- ---- %
  delay = 0.005e-3;
  phase_tx = 0;
  phase_rx = 2*pi*0.13;
  fc_tx = fc;
  fc_rx = fc*1.01;
  SNR = 0;
  % ---- ---- Good Paramters ---- ---- %
  delay = 0.0;
  phase_tx = 0;
  phase_rx = 0;
  fc_tx = fc;
  fc_rx = fc;
  SNR = 100;
  % ---- ---- Pulse shaping ---- ---- %
  pulse_shaping = sin(2*pi*(0:dt:Ts-dt)/(Ts-dt));
  pulse_shaping = ones(1,lt);
% ---- ---- END Transceiver Paramters ---- ---- %



if PLOT_TX
  figure(1);
  hold on;
endif;
if PLOT_RX
  figure(2);
  hold on;
endif;


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
  xi_k = i_seq(mod(k,4)+1);
  xq_k = q_seq(mod(k,4)+1);

  xsi_t = pulse_shaping .* xi_k;
  xsq_t = pulse_shaping .* xq_k;

  xi_t = xsi_t .* cos(2*pi*fc_tx*t+phase_tx);
  xq_t = xsq_t .* sin(2*pi*fc_tx*t+phase_tx);

  x_t = xi_t + xq_t;
  x2_t = [x2_t(1+lt:2*lt) x_t];
  % ---- ---- END Transmitter ---- ---- %




  % ---- ---- Receiver ---- ---- %;
#   y_t = awgn(x_t, SNR);
  delay_idx = mod(floor(delay/dt),lt)+1;
  idx = [delay_idx:(delay_idx+lt-1)];
  y_t = awgn(x2_t(idx), SNR);
  yi_t = y_t .* cos(2*pi*fc_rx*t+phase_rx);
  yq_t = y_t .* sin(2*pi*fc_rx*t+phase_rx);

  % TODO The LBP filter is missing
#   lbp_length = 4;
#   lbp_filter = [ones(1,lbp_length) zeros(1,lt-lbp_length)];
#   lbp_filter /= sum(lbp_filter);
#   yi_t = filter(lbp_filter, [1 zeros(1,19)], [xi_aux_t yi_non_lbp_t])(lt+1:2*lt);
#   yq_t = filter(lbp_filter, [1 zeros(1,19)], [xi_aux_t yq_non_lbp_t])(lt+1:2*lt);
#   xi_aux_t = yi_non_lbp_t;
#   xq_aux_t = yq_non_lbp_t;
  aux = yi_t + j*yq_t;
  aux = fft(aux) .* [1 1 zeros(1,lt-2)];
  aux = ifft(aux);
  yi_filter_t = real(aux);
  yq_filter_t = imag(aux);
  yi_k = yi_filter_t(1);
  yq_k = yq_filter_t(1);
  % ---- ---- END Receiver ---- ---- %




  % ---- ---- Plot Transmitter Signals ---- ---- %
  if (PLOT_TX == 1)
    figure(1);
    plot([t(1) t(1)+Ts],[0,0],'k-');
    stem(t(1), xi_k,'b','linewidth',3);
    stem(t(1), xq_k,'r','linewidth',2);

    plot([t(1) t(1)+Ts],[-3,-3],'k-');
    plot(t,xsi_t-3,'b:','linewidth',2);
    plot(t,xsq_t-3,'r-','linewidth',1);

    plot([t(1) t(1)+Ts],[-6,-6],'k-');
    plot(t,xi_t-6,'b:','linewidth',1);
    plot(t,xq_t-6,'r-','linewidth',1);

    plot([t(1) t(1)+Ts],[-9,-9],'k-');
    plot(t,x_t-9,'g-','linewidth',1);
    figure(1);
    axis([t(1)-window_size*Ts t(1)+Ts -11 2]); %TODO: implement a plotting buffer with a windown from 5 to 10 Ts
    drawnow ("expose");
  endif;
  % ---- ---- END Plot Transmitter Signals ---- ---- %


  % ---- ---- Plot Receiver Signals ---- ---- %
  if (PLOT_RX == 1)
    figure(2);
    plot(t,y_t-9,'g-','linewidth',1)
    plot([t(1) t(1)+Ts],[-9,-9],'k-');
    axis([t(1)-window_size*Ts t(1)+Ts -11 2]); %TODO: implement a plotting buffer with a windown from 5 to 10 Ts
    drawnow ("expose");

    plot([t(1) t(1)+Ts],[-6,-6],'k-');
    plot(t,yi_t-6,'b:','linewidth',1)
    plot(t,yq_t-6,'r-','linewidth',1)

    plot([t(1) t(1)+Ts],[-3,-3],'k-');
    plot(t,yi_filter_t-3,'b:','linewidth',2)
    plot(t,yq_filter_t-3,'r-','linewidth',1)

    plot([t(1) t(1)+Ts],[0,0],'k-');
    stem(t(1), yi_k,'b','linewidth',3);
    stem(t(1), yq_k,'r','linewidth',2);
  endif;
  % ---- ---- END Plot Receiver Signals ---- ---- %


  % ---- ---- Update time and plots ---- ---- %
  t = t + Ts;
  k++;
  fflush(stdout);
  % ---- ---- END Update time and plots ---- ---- %

endwhile;
% ---- ---- END Main Loop ---- ---- %





