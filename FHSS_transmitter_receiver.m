%***************Source*********************
nBits = 20;%number of source bits to transmit
Rb = 20e3; %bit rate of source information in bps
%**********BFSK definitions*****************
fsk_type = 'NONCOHERENT'; %BFSK generation type at the transmitter
h = 1; %modulation index (0.5=coherent BFSK/1= non-coherent BFSK)
Fc = 50e3; %center frequency of BFSK
%**********Frequency Allocation*************
Fbase = 200e3; %The center frequency of the first channel
Fspace = 100e3;%freq. separation between adjacent hopping channels
Fs = 10e6; %sufficiently high sampling frequency for discretization
%*********Frequency Hopper definition*******
G = [1 0 0 1 1]; X=[0 0 0 1];%LFSR generator poly and seed
hopType = 'FAST_HOP'; %FAST_HOP or SLOW_HOP for frequency hopping
%--------Derived Parameters-------
Tb = 1/Rb ;%bit duration of each information bit.
L = Tb*Fs;%num of discrete-time samples in each bit
Fd = h/Tb; %frequency separation of BFSK frequencies
%Adjust num of samples in a hop duration based on hop type selected
if strcmp(hopType,'FAST_HOP'),%hop duration less than bit duration
Lh = L/4; %4 hops in a bit period
nHops = 4*nBits; %total number of Hops during the transmission
else%default set to SLOW_HOP: hop duration more than bit duration
Lh = L*4; %4 bit worth of samples in a hop period
nHops = nBits/4; %total number of Hops during the transmission
end
%-----Simulate the individual blocks----------------------
d = rand(1,nBits) > 0.5 ; %random information bits
[s_m,t,phase,dt]=bfsk_mod(d,Fc,Fd,L,Fs,fsk_type);%BFSK modulation
c = hopping_chip_waveform(G,X,nHops,Lh,Fbase,Fspace,Fs);%Hopping wfm
s = s_m.*c.';%mix BFSK waveform with hopping frequency waveform
n = 0;%Left to the reader -modify for AWGN noise(see prev chapters)
r = s+n; %received signal with noise
v = r.*c.'; %mix received signal with synchronized hopping freq wave
d_cap = bfsk_noncoherent_demod(v,Fc,Fd,L,Fs); %BFSK demod
bie = sum(d ~= d_cap); %Calculate bits in error
disp(['Bits in Error: ', num2str(bie)]);
%--------Plot waveforms at various stages of tx/rx-------
figure;
subplot(2,1,1);plot(t,dt) ;title('Source bits - d(t)');
subplot(2,1,2);plot(t,s_m);title('BFSK modulated - s_m(t)')
figure;
subplot(3,1,1);plot(t,s_m);title('BFSK modulated - s_m(t)')
subplot(3,1,2);plot(t,c);title('Hopping waveform at Tx - c(t)')
subplot(3,1,3);plot(t,s);title('FHSS signal - s(t)')
figure;
subplot(3,1,1);plot(t,r);title('Received signal - r(t)')
subplot(3,1,2);plot(t,c);title('Synced hopping waveform at Rx-c(t)')
subplot(3,1,3);plot(t,v);title('Signal after mixing with hop pattern-v(t)')

function [s,t,phase,at] = bfsk_mod(a,Fc,Fd,L,Fs,fsk_type)
%Function to modulate an incoming binary stream using BFSK
%a - input binary data stream (0's and 1's) to modulate
%Fc - center frequency of the carrier in Hertz
%Fd - frequency separation measured from Fc
%L - number of samples in 1-bit period
%Fs - Sampling frequency for discrete-time simulation
%fsk_type - 'COHERENT' (default) or 'NONCOHERENT' FSK generation
%at each bit period when generating the carriers
%s - BFSK modulated signal
%t - generated time base for the modulated signal
%phase - initial phase generated by modulator, applicable only for
%coherent FSK. It can be used when using coherent detection at Rx
%at - data waveform for the input data
phase=0;
at = kron(a,ones(1,L)); %data to waveform
t = (0:1:length(at)-1)/Fs; %time base
if strcmpi(fsk_type,'NONCOHERENT'),
c1 = cos(2*pi*(Fc+Fd/2)*t+2*pi*rand);%carrier 1 with random phase
c2 = cos(2*pi*(Fc-Fd/2)*t+2*pi*rand);%carrier 2 with random phase
else
phase=2*pi*rand;%random phase from uniform distribution [0,2pi)
c1 = cos(2*pi*(Fc+Fd/2)*t+phase);%carrier 1 with random phase
c2 = cos(2*pi*(Fc-Fd/2)*t+phase);%carrier 2 with random phase
end
s = at.*c1 +(-at+1).*c2; %BFSK signal (MUX selection)
doPlot=0;
if doPlot,
figure;subplot(2,1,1);plot(t,at);subplot(2,1,2);plot(t,s);
end;
end


function a_cap = bfsk_noncoherent_demod(r,Fc,Fd,L,Fs)
%Non-coherent demodulation of BFSK modulated signal
%r - BFSK modulated signal at the receiver
%Fc - center frequency of the carrier in Hertz
%Fd - frequency separation measured from Fc
%L - number of samples in 1-bit period
%Fs - Sampling frequency for discrete-time simulation
%a_cap - data bits after demodulation
t = (0:1:length(r)-1)/Fs; %time base
F1 = (Fc+Fd/2); F2 = (Fc-Fd/2);
%define four basis functions
p1c = cos(2*pi*F1*t); p2c = cos(2*pi*F2*t);
p1s = -sin(2*pi*F1*t); p2s = -sin(2*pi*F2*t);
%multiply and integrate from 0 to Tb
r1c = conv(r.*p1c,ones(1,L)); r2c = conv(r.*p2c,ones(1,L));
r1s = conv(r.*p1s,ones(1,L)); r2s = conv(r.*p2s,ones(1,L));
%sample at every sampling instant
r1c = r1c(L:L:end); r2c = r2c(L:L:end);
r1s = r1s(L:L:end); r2s = r2s(L:L:end);
x = r1c.^2 + r1s.^2; y = r2c.^2 + r2s.^2;%square and add
a_cap=(x-y)>0; %compare and decide
end


function freqTable = gen_FH_code_table(Fbase,Fspace,Fs,Lh,N)
%Generate frequency translation table for Frequency Hopping (FH)
%Fbase - base frequency (Hz) of the hop
%Fspace - channel spacing (Hz) of the hop
%Fs - sampling frequency (Hz)
%Lh - num of discrete time samples in each hop period
%N - num of total hopping frequencies required(full period of LFSR)
%Return the frequency translation table
t = (0:1:Lh-1)/Fs; %time base for each hopping period
freqTable = zeros(N,Lh);%Table to store N different freq waveforms
for i=1:N, %generate frequency translation table for all N states
Fi=Fbase+(i-1)*Fspace;
freqTable(i,:) = cos(2*pi*Fi*t);
end
end


function [c] = hopping_chip_waveform(G,X,nHops,Lh,Fbase,Fspace,Fs)
%Generate Frequency Hopping chip sequence based on m-sequence LFSR
%G,X - Generator poly. and initial seed for underlying m-sequence
%nHops - total number of hops needed
%Lh - number of discrete time samples in each hop period
%Fbase - base frequency (Hz) of the hop
%Fspace - channel spacing (Hz) of the hop
%Fs - sampling frequency (Hz)
[prbs,STATES] = lfsr( G, X); %initialize LFSR
N = length(prbs); %PRBS period
LFSRStates= repeatSequence(STATES.',nHops);%repeat LFSR states depending on nHops
freqTable=gen_FH_code_table(Fbase,Fspace,Fs,Lh,N); %freq translation
c = zeros(nHops,Lh); %place holder for the hopping sequence waveform
for i=1:nHops,%for each hop choose one freq wave based on LFSR state
LFSRState = LFSRStates(i); %given LFSR state
c(i,:) = freqTable(LFSRState,:);%choose corresponding freq. wave
end
c=c.';c=c(:); %transpose and serialize as a single waveform vector
end


function [y,states] = lfsr(G,X)
%Galois LFSR for m-sequence generation
% G - polynomial (primitive for m-sequence) arranged as vector
% X - initial state of the delay elements arranged as vector
% y - output of LFSR
% states - gives the states through which the LFSR has sequenced.
%This is particularly helpful in frequency hopping applications
%The function outputs the m-sequence for a single period
%Sample call:
%3rd order LFSR polynomial:xˆ5+xˆ2+1=>g5=1,g4=0,g3=0,g2=1,g1=0,g0=1
%with intial states [0 0 0 0 1]: lfsr([1 0 0 1 0 1],[0 0 0 0 1])
g=G(:); x=X(:); %serialize G and X as column vectors
if length(g)-length(x)~=1,
error('Length of initial seed X0 should be equal to the number of delay elements (length(g)-1)');
end
%LFSR state-transistion matrix construction
L = length(g)-1; %order of polynomial
A0 = [zeros(1,L-1); eye(L-1)]; %A-matrix construction
g=g(1:end-1);
A = [A0 g]; %LFSR state-transistion matrix
N = 2^L-1; %period of maximal length sequence
y = zeros(1,length(N));%array to store output
states=zeros(1,length(N));%LFSR states(useful for Frequeny Hopping)
for i=1:N, %repeate for each clock period
states(i)=bin2dec(char(x.'+'0'));%convert LFSR states to a number
y(i) = x(end); %output the last bit
x = mod(A*x,2); %LFSR equation
end
end

function [y]= repeatSequence(x,N)
%Repeat a given sequence x of arbitrary length to match the given
%length N. This function is useful to repeat any sequence
%(say PRBS) to match the length of another sequence
x=x(:); %serialize
xLen = length(x); %length of sequence x
%truncate or extend sequence x to suite the given length N
if xLen >= N, %Truncate x when sequencelength less than N
y = x(1:N);
else
temp = repmat(x,fix(N/xLen),1); %repeat sequence integer times
residue = mod(N,xLen); %reminder when dividing N by xLen
%append reminder times
if residue ~=0, temp = [temp; x(1:residue)]; end
y = temp; %repeating sequence matching length N
end
end

