% Define parameters
Rb = 100; Rc = 1000; L = 32;
prbsType = 'GOLD';
G1 = [1 1 1 1 0 0 0 1]; G2 = [1 0 0 1 0 0 0 1];
X1 = [0 0 0 0 0 0 1]; X2 = [0 0 0 0 0 0 1];
d = rand(1,2) >= 0.5;

% Transmitter - Generate DSSS signal
[s_t, carrier_ref, prbs_ref] = dsss_transmitter(d, prbsType, G1, G2, X1, X2, Rb, Rc, L);

% Plot Transmitter Output in Frequency Domain
S_f = abs(fftshift(fft(s_t))); % Fourier Transform of Transmitter Output
f = linspace(-0.5, 0.5, length(S_f)); % Frequency vector for plotting

figure;
subplot(2,1,1);
plot(f, S_f);
title('Transmitter Output in Frequency Domain');
xlabel('Normalized Frequency');
ylabel('Magnitude');

% Receiver - Simulate received signal and demodulate
r_t = s_t; % Assume ideal channel (no noise or fading)
d_cap = dsss_receiver(r_t, carrier_ref, prbs_ref, Rb, Rc, L);

% Plot Receiver Output in Frequency Domain
V_f = abs(fftshift(fft(r_t))); % Fourier Transform of Receiver Output
subplot(2,1,2);
plot(f, V_f);
title('Receiver Output in Frequency Domain');
xlabel('Normalized Frequency');
ylabel('Magnitude');

function [prbs] = generatePRBS(prbsType,G1,G2,X1,X2)
%Generate PRBS sequence - choose from either msequence or gold code
%prbsType - type of PRBS generator - 'MSEQUENCE' or 'GOLD'
% If prbsType == 'MSEQUENCE' G1 is the generator poly for LFSR
% and X1 its seed. G2 and X2 are not used
% If prbsType == 'GOLD' G1,G2 are the generator polynomials
% for LFSR1/LFSR2 and X1,X2 are their initial seeds.
%G1,G2 - Generator polynomials for PRBS generation
%X1,X2 - Initial states of LFSRs
%The PRBS generators results in 1 period PRBS,
%need to repeat it to suit the data length
if strcmpi(prbsType,'MSEQUENCE'),
prbs=lfsr( G1, X1);%use only one poly and initial state vector
elseif strcmpi(prbsType,'GOLD'),
prbs=gold_code_generator(G1,G2,X1,X2);%full length Gold sequence
else %Gold codes as default
G1=[1 1 1 1 0 1]; G2 = [1 0 0 1 0 1]; %LFSR polynomials
X1 = [ 0 0 0 0 1]; X2=[ 0 0 0 0 1] ; %initial state of LFSR
prbs = gold_code_generator(G1,G2,X1,X2);
end
end

function [y] = gold_code_generator( G1,G2,X1,X2)
%Implementation of Gold code generator
%G1-preferred polynomial 1 for LFSR1 arranged as [g0 g1 g2 ... gL-1]
%G2-preferred polynomial 2 for LFSR2 arranged as [g0 g1 g2 ... gL-1]
%X1-initial seed for LFSR1 [x0 x1 x2 ... xL-1]
%X2-initial seed for LFSR2 [x0 x1 x2 ... xL-1]
%y-Gold code
%The function outputs the m-sequence for a single period
%Sample call:
% 7th order preferred polynomials [1,2,3,7] and [3,7] (polynomials : 1+x+xˆ2+xˆ3+x7 and 1+xˆ3+xˆ7)
% i.e, G1 = [1,1,1,1,0,0,0,1] and G2=[1,0,0,1,0,0,0,1]
% with intial states X1=[0,0,0,0,0,0,0,1], X2=[0,0,0,0,0,0,0,1]:
%gold_code_generator(G1,G2,X1,X2)
g1=G1(:); x1=X1(:); %serialize G1 and X1 matrices
g2=G2(:); x2=X2(:); %serialize G2 and X2 matrices
if length(g1)~=length(g2) && length(x1)~=length(x2),
    error('Length mismatch between G1 & G2 or X1 & X2');
end
%LFSR state-transistion matrix construction
L = length(g1)-1; %order of polynomial
A0 = [zeros(1,L-1); eye(L-1)]; %A-matrix construction
g1=g1(1:end-1);g2=g2(1:end-1);
A1 = [A0 g1]; %LFSR1 state-transistion matrix
A2 = [A0 g2]; %LFSR2 state-transistion matrix
N = 2^L-1; %period of maximal length sequence
y = zeros(1,length(N));%array to store output
for i=1:N, %repeate for each clock period
y(i)= mod(x1(end)+x2(end),2);%XOR of outputs of LFSR1 & LFSR2
x1 = mod(A1*x1,2); %LFSR equation
x2 = mod(A2*x2,2); %LFSR equation
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

function [s_t,carrier_ref,prbs_ref]=dsss_transmitter(d,prbsType,G1,G2,X1,X2,Rb,Rc,L)
% Direct Sequence Spread Spectrum (DSSS) transmitter - returns the DSSS
% waveform (s), reference carrier, the prbs reference waveform for use
% in synchronization in the receiver
%d - input binary data stream
%prbsType - type of PRBS generator - 'MSEQUENCE' or 'GOLD'
% If prbsType == 'MSEQUENCE' G1 is the generator poly for LFSR
% and X1 its seed. G2 and X2 are not used
% If prbsType == 'GOLD' G1,G2 are the generator polynomials
% for LFSR1/LFSR2 and X1,X2 are their initial seeds.
%G1,G2 - Generator polynomials for PRBS generation
%X1,X2 - Initial states of LFSRs
%Rb - data rate (bps) for the data d
%Rc - chip-rate (Rc >> Rb AND Rc is integral multiple of Rb)
%L - oversampling factor for waveform generation
prbs = generatePRBS(prbsType,G1,G2,X1,X2);
prbs=prbs(:); d=d(:); %serialize
dataLen= length(d)*(Rc/Rb);%required PRBS length to cover the data
prbs_ref= repeatSequence(prbs,dataLen);%repeat PRBS to match data
d_t = kron(d,ones(L*Rc/Rb,1)); %data waveform
prbs_t = kron(prbs_ref,ones(L,1)); %spreading sequence waveform
sbb_t = 2*xor(d_t,prbs_t)-1; %XOR data and PRBS, convert to bipolar
n=(0:1:length(sbb_t)-1).'; carrier_ref=cos(2*pi*2*n/L);
s_t = sbb_t.*carrier_ref; %modulation,2 cycles per chip
figure(1); %Plot waveforms
subplot(3,1,1); plot(d_t); title('data sequence'); hold on;
subplot(3,1,2); plot(prbs_t); title('PRBS sequence');
subplot(3,1,3); plot(s_t); title('DS-SS signal (baseband)'); end

function d_cap = dsss_receiver(r_t,carrier_ref,prbs_ref,Rb,Rc,L)
%Direct Sequence Spread Spectrum (DSSS) Receiver (Rx)
%r_t - received DS-SS signal from the transmitter (Tx)
%carrier_ref - reference carrier (synchronized with transmitter)
%prbs_ref - reference PRBS signal(synchronized with transmitter)
%Rb - data rate (bps) for the data d
%Rc - chip-rate ((Rc >> Rb AND Rc is integral multiple of Rb)
%L - versampling factor used for waveform generation at the Tx
%The function first demodulates the receiver signal using the
%reference carrier and then XORs the result with the reference
%PRBS. Finally returns the recovered data.
%------------BPSK Demodulation----------
v_t = r_t.*carrier_ref;
x_t=conv(v_t,ones(1,L)); %integrate for Tc duration
y = x_t(L:L:end);%sample at every Lth sample (i.e, every Tc instant)
z = ( y > 0 ).'; %Hard decision (gives demodulated bits)
%-----------De-Spreading----------------
y = xor(z,prbs_ref.');%reverse the spreading process using PRBS ref.
d_cap = y(Rc/Rb:Rc/Rb:end); %sample at every Rc/Rb th symbol
end

