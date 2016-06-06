clear all;
close all;
%clf;

%---------------- Initial Parameters ----------------

n = str2num (input ("Enter the number of bits in the bitstream\n", "s")); %number of bits in the FIFO
choice = menu ("Which QAM do you want to use", "4", "16", "64");  %number of symbols in the constellation
if choice == 1
  M = 4;
elseif choice == 2
  M = 16;
elseif choice == 3
  M = 64;
endif
k = log2(M); %number of bits in one symbol
w = 2*pi*10;  % angular frequency

%----------------Creation of the bitstream ---------

dataIn = randi(2,n,1) - 1;  % Generate vector of binary data

stem(dataIn,"filled");
title('Random Bits');
xlabel('Bit Index');
ylabel('Binary Value');

%----------------- QAM Mapper --------------------

numberOfSymbols = length(dataIn)/k; %number of symbols in the FIFO

symbol = 0;

for i = 1:n % rewshape dataIn into dataInPack with one case for one symbol composed of k bits 

  if mod(i,k) == 0
    dataInPack(i/k) = symbol*10 + dataIn(i);
    symbol = 0;
  else
    symbol = symbol*10 + dataIn(i);
   endif
   
endfor

%dataInPack

for i = 1:numberOfSymbols  % symbols : matrix which will contain the value of I(first column) and Q (second column) for each symbol

  symbolDe  = bin2dec(mat2str(dataInPack(i)));
  
  symbols(i,1) = 0;
  symbols(i,2) = 0;
  
  if M == 64
  
    if mod((symbolDe+2),8) == 0
      symbols(i,1) = +1;
    elseif mod((symbolDe+1),8) == 0
      symbols(i,1) = +3;
    elseif mod((symbolDe+3),8) == 0
      symbols(i,1) = +5;
    elseif mod((symbolDe+4),8) == 0
      symbols(i,1) = +7;
    elseif mod((symbolDe+6),8) == 0
      symbols(i,1) = -1;
    elseif mod((symbolDe+5),8) == 0
      symbols(i,1) = -3;
    elseif mod((symbolDe+7),8) == 0
      symbols(i,1) = -5;
    elseif mod(symbolDe,8) == 0
      symbols(i,1) = -7;
    endif
  
    if and(16 <= symbolDe, symbolDe < 24)
      symbols(i,2) = +1;
    elseif and(24 <= symbolDe, symbolDe < 32)
      symbols(i,2) = +3;
    elseif and(8 <= symbolDe, symbolDe < 16)
      symbols(i,2) = +5;
    elseif and(0 <= symbolDe, symbolDe < 8)
      symbols(i,2) = +7;
    elseif and(48 <= symbolDe, symbolDe < 56)
      symbols(i,2) = -1;
    elseif and(56 <= symbolDe, symbolDe < 64)
      symbols(i,2) = -3;
    elseif and(40 <= symbolDe, symbolDe < 48)
      symbols(i,2) = -5;
    elseif and(32 <= symbolDe, symbolDe < 40)
      symbols(i,2) = -7;
    endif
      
  elseif M == 16
  
    if mod(symbolDe,4) == 0
      symbols(i,2) = -3;
    elseif mod((symbolDe + 2),4) == 0
      symbols(i,2) = +3; 
    elseif mod((symbolDe + 3),4) == 0
      symbols(i,2) = -1;
    elseif mod((symbolDe + 1),4) == 0
      symbols(i,2) = +1;
    endif
      
    if symbolDe < 4
      symbols(i,1) = -3;
    elseif and(4 <= symbolDe, symbolDe < 8)
      symbols(i,1) = -1;
    elseif and(8 <= symbolDe, symbolDe < 12)
      symbols(i,1) = +3;
    elseif and(12 <= symbolDe, symbolDe < 16)
      symbols(i,1) = +1;
    endif
    
  elseif M == 4
    if mod(symbolDe,2) == 0
      symbols(i,1) = +1;
    elseif mod((symbolDe + 1),2) == 0
      symbols(i,1) = -1;
    endif
      
    if symbolDe < 2
      symbols(i,2) = +1;
    elseif 1 < symbolDe < 4
      symbols(i,2) = -1;
    endif   
  endif

endfor

%symbols

%------------------------------Construction and plot of the signal-----------

signal = symbols(:,1) + j*symbols(:,2);

scatterplot(signal);
title('Constellation before sending');
grid;
xlabel('I');
ylabel('Q');
if M == 4
  axis([-2 2 -2 2]);
elseif M == 16
  axis([-4 4 -4 4]);
elseif M == 64
  axis([-8 8 -8 8]);
endif

t = 0:1/k:length(dataIn)-1/k;
modulation = real(signal)*cos(w*t) + imag(signal)*sin(w*t);

%--------------------------------Channel with noise addition----------------

EbNo = 10;
numSamplesPerSymbol = 1;    % Oversampling factor
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);

signalAfter = awgn(signal,snr,'measured');

%----------------------------- Deconstruction of the signal and plot of the received signal--------------

scatterplot(signalAfter);
title('Constellation after sending');
grid;
xlabel('I');
ylabel('Q');
if M == 4
  axis([-2 2 -2 2]);
elseif M == 16
  axis([-4 4 -4 4]);
elseif M == 64
  axis([-8 8 -8 8]);
endif

phi = 0;

signalAfterDeconstructionI = signalAfter*cos(w*t + phi);
signalAfterDeconstructionQ = signalAfter*sin(w*t + phi);

%-------------------------------FIR-------------------------------------------