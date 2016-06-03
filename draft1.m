clear all;
close all;
clf;

%---------------- Initial Parameters ----------------

n = 10; %number of bits in the FIFO
M = 4;  %number of symbols in the constellation
k = log2(M); %number of bits in one symbol


%----------------Creation of the bitstream ---------

dataIn = randi(2,n,1) - 1;  % Generate vector of binary data

stem(dataIn,"filled");
title('Random Bits');
xlabel('Bit Index');
ylabel('Binary Value');

%------------- QAM Mapper --------------------

dataInPack = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary k-tuples, k = log2(M)
dataSymbolsIn = bi2de(dataInPack);                 % Convert to integers


figure; % Create new figure window.
stem(dataSymbolsIn(1:10));
title('Random Symbols');
xlabel('Symbol Index');
ylabel('Integer Value');
