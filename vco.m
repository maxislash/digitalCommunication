 % ---- ---- Simulation Paramters ---- ---- %
Ts = 1e-3;
dt = 1e-6;
fs = 1/dt;
end_t = 1e-3;
t = 0:dt:Ts-dt;
lt = length(t);
fc = 1e3;
delay = 0.0;
phase_rx = 0;
phase_tx = 0;
fc_tx = fc;
fc_rx = fc;
SNR = 100;   

% first order pll
alpha = 0.01;

% VCO voltage which controlls the frequency. at v=0 it's exactly at f.
v = 0;

% cosine output
c = 1;
% delayed cosine by one timestep
c_delay = 0;
% sine output
s = 0;
% delayed sine by one timestep
s_delay = 0;

% data from the VCO for debugging purposes
sine = zeros(Ts,1);
cosine =zeros(Ts,1); 
vco = zeros(Ts,1); 

for step = 1:dt/dt:Ts/dt
  
  % this is part of the PLL
  % "voltage" controlled oscillator 
  f0 = fc_rx/fs + v*alpha;
  c_delay = c;
  s_delay = s;
  c = c_delay * cos(2*pi*f0) - s_delay * sin(2*pi*f0);
  s = s_delay * cos(2*pi*f0) + c_delay * sin(2*pi*f0);
  % let's save everything in handy vectors for plotting
  sine(step) = s;
  cosine(step) = c;
  %qvco(step) = v;
  % end VCO

  % save our carriers
  % !!! 90 degree phase shift so the sine becomes the inphase
  %     signal and the cosine the quadrature signal
  carrier_inph(step) = -s;
  carrier_quad(step) = c;

endfor

figure(1);

subplot(2,1,1);
plot(carrier_inph,'LineWidth',2);
xlabel('Samples');
ylabel('Amplitude');
title('inph');
grid on;

subplot(2,1,2);
plot(carrier_quad,'LineWidth',2);
xlabel('Samples');
ylabel('Amplitude');
title('quad');
grid on;
