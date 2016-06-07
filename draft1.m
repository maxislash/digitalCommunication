clear all;
close all;
%clf;

%---------------- Initial Parameters ----------------

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
k = log2(M); %number of bits in one symbol
w = 2*pi*10;  % angular frequency

if mod(n, k) ~= 0
    error('numberOfBits must be a multiple of log2(M).');
end

%----------------Creation of the bitstream ---------

dataIn = randi(2,n,1) - 1;  % Generate vector of binary data

%stem(dataIn,"filled");
%title('Random Bits');
%xlabel('Bit Index');
%ylabel('Binary Value');

%----------------- QAM Mapper --------------------

numberOfSymbols = length(dataIn)/k; %number of symbols in the FIFO

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

mappedSymbols = zeros(1, numberOfSymbols);

% Map bits to symbols
for i = 1:k:length(dataIn)
   
    symbolBits = dataIn(i:i+k-1);
    
    if M == 4
      symbolIndex = 2^1 * symbolBits(1) + 2^0 * symbolBits(2);
    elseif M == 16
      symbolIndex = 2^3 * symbolBits(1) + 2^2 * symbolBits(2) + 2^1 * symbolBits(3) + 2^0 * symbolBits(4);
    elseif M == 64
      symbolIndex = 2^5 * symbolBits(1) + 2^4 * symbolBits(2) + 2^3 * symbolBits(3) + 2^2 * symbolBits(4) +  2^1 * symbolBits(5) + 2^0 * symbolBits(6);
    endif

     % Mapping
    mappedSymbols((i - 1)/k + 1) = mappingTable( symbolIndex + 1);
endfor

%------------------------------Construction and plot of the signal-----------

signal = mappedSymbols;

t = 0:1/k:length(dataIn)/(k*k)-1/k;
modulation = real(signal).*cos(w.*t) + imag(signal).*sin(w.*t);

%--------------------------------Channel with noise addition----------------

EbNo = 10;
numSamplesPerSymbol = 1;    % Oversampling factor
snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);

signalAfter = awgn(signal,snr,'measured');

%----------------------------- Deconstruction of the signal and plot of the received signal--------------

sPlotFig = scatterplot(signal,1,0,'k*');
hold on
scatterplot(signalAfter,1,0,'g.', sPlotFig);
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
title('Signal before and after sending');
legend('Before Sending','After Sending');

phi = 0;

signalAfterDeconstructionI = signalAfter.*cos(w*t + phi);
signalAfterDeconstructionQ = signalAfter.*sin(w*t + phi);

%-------------------------------FIR-------------------------------------------

% Low pass filtering with a Butterworth filter
[b,a]=butter(2,0.4);
Yi=filter(b,a,signalAfterDeconstructionI);
Yk=filter(b,a,signalAfterDeconstructionQ);

%----------------------------- QAM demapper -------------------------------

%receivedSignal = Yi + j*Yk;
%for i = 1:length(receivedSignal)
%  [mindiff minIndex] = min(receivedSignal(i) - mappingTable);
%  symbolIndex = minIndex - 1;
%  bitString = dec2bin(symbolIndex, 4);
%end
%
%  receivedBits((i-1)*4 + 1) = str2double(bitString(1));
%  receivedBits((i-1)*4 + 2) = str2double(bitString(2));
%  receivedBits((i-1)*4 + 3) = str2double(bitString(3));
%  receivedBits((i-1)*4 + 4) = str2double(bitString(4));
%receivedBits = zeros(1, numberOfBits / 4);

