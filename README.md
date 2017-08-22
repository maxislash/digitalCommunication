# digitalCommunication

Designed a digital communication system using QAM modulation on Octave with Symbol Timing and Carrier Frequency Recovery. 

QAM is useful in digital communication because the quadrature carriers create a combination of phase-shift and amplitude-shift keying which is less subject to noise. 
The transmitter is composed of a QAM Mapper, pulse shaping, modulator and a sum block. 
The receiver is composed of a demodulator, low pass filter, symbol timing recovery block carrier frequency and phase recovery block and QAM demapper. 
The Symbol Timing Recovery was made using an Early/Late Symbol Recovery algorithm and the frequency and phase recovery used a Costas loop. 

The system is coded in Matlab/Octave.
