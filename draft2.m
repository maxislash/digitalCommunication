close all;
clear all;
graphics_toolkit ("gnuplot");

rx_p = 0;
rx_f = 1e3;

vco_p(1) = 4;
vco_f = rx_f*1.01321;

delta_t = 1e-5;
max_t = 5e-1;
t = 0:delta_t:max_t;

rx_t = cos(2*pi*rx_f*t + rx_p);


maf_l = 500;
fifo = zeros(1,maf_l);

maf_f_l = 100;
fifo_f = zeros(1,maf_f_l);

maf(1) = 0;
maf_f(1) = 0;

for i = 2:length(rx_t)

  vco_t(i) = -sin(2*pi*vco_f(i-1)*t(i) + vco_p(i-1));
  aux_t(i) = rx_t(i) * vco_t(i);
  maf(i) = maf(i-1) + aux_t(i);
  maf(i) -= fifo(maf_l);

  fifo(2:maf_l) = fifo(1:maf_l-1);
  fifo(1) = aux_t(i);

  vco_p(i) = vco_p(i-1) + maf(i-1)/maf_l/100;

  maf_f(i) = maf_f(i-1) + (vco_p(i) - vco_p(i-1))/delta_t;
  maf_f(i) -= fifo_f(maf_f_l);
  fifo_f(2:maf_f_l) = fifo_f(1:maf_f_l-1);
  fifo_f(1) = (vco_p(i) - vco_p(i-1))/delta_t;

  if (mod(i,10000) == 0)
    vco_f(i) = vco_f(i-1) + maf_f(i)/pi/2/maf_f_l;
  else
    vco_f(i) = vco_f(i-1);
  endif;

endfor;

p_step = 100;

figure(1); hold on;
plot(t(1:p_step:end),aux_t(1:p_step:end));
plot(t(1:p_step:end),maf(1:p_step:end)/maf_l*2,'r;phase error;');
plot(t(1:p_step:end),vco_p(1:p_step:end),'g;vco phase;');
# plot(t(1:p_step:end),vco_f(1:p_step:end)-rx_f(1:p_step:end),'c;vco f - rx f;')
plot(t(1:p_step:end),maf_f(1:p_step:end)/1000/pi,'k;freq error;')

figure(3); hold on;
plot(t(1:p_step:end),rx_f*ones(1,length(t(1:p_step:end))),'b;rx f;');
plot(t(1:p_step:end),vco_f(1:p_step:end),'c;vco f;');
# plot(t,maf_f/10000,'k;freq error;');

# figure(2); hold on;
# plot(t,rx_t);
# plot(t,vco_t,'r');


