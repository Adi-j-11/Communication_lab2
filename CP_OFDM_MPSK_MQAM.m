clearvars; clc;
%--------Simulation parameters----------------
nSym = 10^4; % Number of OFDM Symbols to transmit
EbN0dB = 0:2:20; % Bit to noise ratio
% Define modulation types and orders
PSK_MOD_TYPES = {'MPSK', 'MPSK', 'MPSK', 'MPSK', 'MPSK', 'MPSK'};
PSK_Ms = [2, 4, 8, 16, 32, 64]; % PSK Modulation orders
QAM_MOD_TYPES = {'MQAM', 'MQAM', 'MQAM', 'MQAM'};
QAM_Ms = [4, 16, 64, 256]; % QAM Modulation orders
N = 64; % FFT size or total number of subcarriers (used + unused)
Ncp = 16; % Number of symbols in the cyclic prefix

% Prepare plots
figure;

% Plot PSK
subplot(2, 1, 1);
hold on;
for idx = 1:length(PSK_Ms)
    M = PSK_Ms(idx);
    MOD_TYPE = PSK_MOD_TYPES{idx};
    
    %--------Derived Parameters--------------------
    k = log2(M); % Number of bits per modulated symbol
    EsN0dB = 10*log10(k) + EbN0dB; % Convert to symbol energy to noise ratio
    errors = zeros(1, length(EbN0dB)); % To store symbol errors
    
    for i = 1:length(EbN0dB) % Monte Carlo Simulation
        for j = 1:nSym
            %-----------------Transmitter--------------------
            d = ceil(M .* rand(1, N)); % Uniformly distributed random symbols from 1:M
            [X, ref] = modulation_mapper(MOD_TYPE, M, d);
            x = ifft(X, N); % IDFT
            s = add_cyclic_prefix(x, Ncp); % Add CP
            
            %-------------- Channel ----------------
            r = add_awgn_noise(s, EsN0dB(i)); % Add AWGN noise r = s + n
            
            %-----------------Receiver----------------------
            y = remove_cyclic_prefix(r, Ncp, N); % Remove CP
            Y = fft(y, N); % DFT
            [~, dcap] = iqOptDetector(Y, ref); % Demapper using IQ detector
            
            %----------------Error counter------------------
            numErrors = sum(d ~= dcap); % Count number of symbol errors
            errors(i) = errors(i) + numErrors; % Accumulate symbol errors
        end
    end
    
    simulatedSER = errors / (nSym * N);
    
    % Plot simulated SER for PSK
    plot(EbN0dB, log10(simulatedSER), 'o-', 'DisplayName', sprintf('%d-PSK', M));
end
title('Simulated SER for PSK Modulation Schemes');
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate (log10)');
legend('show');
grid on;

% Plot QAM
subplot(2, 1, 2);
hold on;
for idx = 1:length(QAM_Ms)
    M = QAM_Ms(idx);
    MOD_TYPE = QAM_MOD_TYPES{idx};
    
    %--------Derived Parameters--------------------
    k = log2(M); % Number of bits per modulated symbol
    EsN0dB = 10*log10(k) + EbN0dB; % Convert to symbol energy to noise ratio
    errors = zeros(1, length(EbN0dB)); % To store symbol errors
    
    for i = 1:length(EbN0dB) % Monte Carlo Simulation
        for j = 1:nSym
            %-----------------Transmitter--------------------
            d = ceil(M .* rand(1, N)); % Uniformly distributed random symbols from 1:M
            [X, ref] = modulation_mapper(MOD_TYPE, M, d);
            x = ifft(X, N); % IDFT
            s = add_cyclic_prefix(x, Ncp); % Add CP
            
            %-------------- Channel ----------------
            r = add_awgn_noise(s, EsN0dB(i)); % Add AWGN noise r = s + n
            
            %-----------------Receiver----------------------
            y = remove_cyclic_prefix(r, Ncp, N); % Remove CP
            Y = fft(y, N); % DFT
            [~, dcap] = iqOptDetector(Y, ref); % Demapper using IQ detector
            
            %----------------Error counter------------------
            numErrors = sum(d ~= dcap); % Count number of symbol errors
            errors(i) = errors(i) + numErrors; % Accumulate symbol errors
        end
    end
    
    simulatedSER = errors / (nSym * N);
    
    % Plot simulated SER for QAM
    plot(EbN0dB, log10(simulatedSER), 's-', 'DisplayName', sprintf('%d-QAM', M));
end
title('Simulated SER for QAM Modulation Schemes');
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate (log10)');
legend('show');
grid on;


function [X,ref]=modulation_mapper(MOD_TYPE,M,d)
%Modulation mapper for OFDM transmitter
% MOD_TYPE - 'MPSK' or 'MQAM' modulation
% M - modulation order, For BPSK M=2, QPSK M=4, 256-QAM M=256 etc..,
% d - data symbols to be modulated drawn from the set {1,2,...,M}
%returns
% X - modulated symbols
% ref -ideal constellation points that could be used by IQ detector
if strcmpi(MOD_TYPE,'MPSK'),
[X,ref]=mpsk_modulator(M,d);%MPSK modulation
else
if strcmpi(MOD_TYPE,'MQAM'),
[X,ref]=mqam_modulator(M,d);%MQAM modulation
else
error('Invalid Modulation specified');
end
end;end

function [s,ref]=mqam_modulator(M,d)
%Function to MQAM modulate the vector of data symbols - d
%[s,ref]=mqam_modulator(M,d) modulates the symbols defined by the vector d
% using MQAM modulation, where M specifies order of M-QAM modulation and
% vector d contains symbols whose values range 1:M. The output s is modulated
% output and ref represents reference constellation that can be used in demod
if(((M~=1) && ~mod(floor(log2(M)),2))==0), %M not a even power of 2
error('Only Square MQAM supported. M must be even power of 2');
end
ref=constructQAM(M); %construct reference constellation
s=ref(d); %map information symbols to modulated symbols
end

function [ref,varargout]= constructQAM(M)
%Function to construct gray codes symbol constellation for M-QAM
% [ref]=constructQAM(M) - returns the ideal signaling points (ref) in a
% symmetric rectangular M-QAM constellation, where M is the level of QAM
% modulation. The returned constellation points are arranged such that the
% index of the points are arranged in a Gray-coded manner. When plotted,
% indices of constellation points will differ by 1 bit.
%
% [ref,I,Q]=constructQAM(M) - returns the ideal signaling points (ref) along
% with the IQ components breakup in a symmetric rectangular M-QAM constellation,
% where M is the level of QAM modulation. The returned constellation points are
% arranged such that the index of the points are arranged in a Gray-coded manner.
n=0:1:M-1; %Sequential address from 0 to M-1 (1xM dimension)
%------Addresses in Kmap - Gray code walk---------------
a=dec2gray(n); %Convert linear addresses to gray code
N=sqrt(M); %Dimension of K-Map - N x N matrix
a=reshape(a,N,N).'; %NxN gray coded matrix
evenRows=2:2:size(a,1); %identify alternate rows
a(evenRows,:)=fliplr(a(evenRows,:));%Flip rows - KMap representation
nGray=reshape(a.',1,M); %reshape to 1xM - Gray code walk on KMap
%Construction of ideal M-QAM constellation from sqrt(M)-PAM
D=sqrt(M); %Dimension of PAM constellation
x=floor(nGray/D);
y=mod(nGray,D);
Ax=2*(x+1)-1-D; %PAM Amplitudes 2m-1-D - real axis
Ay=2*(y+1)-1-D; %PAM Amplitudes 2m-1-D - imag axis
ref=Ax+1i*Ay; %assign complex numbers to reference
if nargout==2 %if only one variable argument is given
varargout{1}=Ax; %Real part (I)
elseif nargout==3 %if two variable arguments are given
varargout{1}=Ax; %Real part (I)
varargout{2}=Ay; %Imaginary part (Q)
end
end

function [s,ref]=mpsk_modulator(M,d)
%Function to MPSK modulate the vector of data symbols - d
%[s,ref]=mpsk_modulator(M,d) modulates the symbols defined by the
%vector d using MPSK modulation, where M specifies the order of
%M-PSK modulation and the vector d contains symbols whose values
%in the range 1:M. The output s is the modulated output and ref
%represents the reference constellation that can be used in demod
ref_i= 1/sqrt(2)*cos(((1:1:M)-1)/M*2*pi);
ref_q= 1/sqrt(2)*sin(((1:1:M)-1)/M*2*pi);
ref = ref_i+1i*ref_q;
s = ref(d); %M-PSK Mapping
end

function s = add_cyclic_prefix(x,Ncp)
%function to add cyclic prefix to the generated OFDM symbol x that
%is generated at the output of the IDFT block
% x - ofdm symbol without CP (output of IDFT block)
% Ncp-num. of samples at x's end that will copied to its beginning
% s - returns the cyclic prefixed OFDM symbol
s = [x(end-Ncp+1:end) x]; %Cyclic prefixed OFDM symbol
end

function [r,n,N0] = add_awgn_noise(s,SNRdB,L)
%Function to add AWGN to the given signal
%[r,n,N0]= add_awgn_noise(s,SNRdB) adds AWGN noise vector to signal
%'s' to generate a %resulting signal vector 'r' of specified SNR
%in dB. It also returns the noise vector 'n' that is added to the
%signal 's' and the spectral density N0 of noise added
%
%[r,n,N0]= add_awgn_noise(s,SNRdB,L) adds AWGN noise vector to
%signal 's' to generate a resulting signal vector 'r' of specified
%SNR in dB. The parameter 'L' specifies the oversampling ratio used
%in the system (for waveform simulation). It also returns the noise
%vector 'n' that is added to the signal 's' and the spectral
%density N0 of noise added
s_temp=s;
if iscolumn(s), s=s.'; end; %to return the result in same dim as 's'
gamma = 10^(SNRdB/10); %SNR to linear scale
if nargin==2, L=1; end %if third argument is not given, set it to 1
if isvector(s),
P=L*sum(abs(s).^2)/length(s);%Actual power in the vector
else %for multi-dimensional signals like MFSK
P=L*sum(sum(abs(s).^2))/length(s); %if s is a matrix [MxN]
end
N0=P/gamma; %Find the noise spectral density
if(isreal(s)),
n = sqrt(N0/2)*randn(size(s));%computed noise
else
n = sqrt(N0/2)*(randn(size(s))+1i*randn(size(s)));%computed noise
end
r = s + n; %received signal
if iscolumn(s_temp), r=r.'; end;%return r in original format as s
end

function y = remove_cyclic_prefix(r,Ncp,N)
%function to remove cyclic prefix from the received OFDM symbol r
% r - received ofdm symbol with CP
% Ncp - num. of samples at beginning of r that need to be removed
% N - number of samples in a single OFDM symbol
% y - returns the OFDM symbol without cyclic prefix
y=r(Ncp+1:N+Ncp);%cut from index Ncp+1 to N+Ncp
end

function [idealPoints,indices]= iqOptDetector(received,ref)
%Optimum Detector for 2-dim. signals (MQAM,MPSK,MPAM) in IQ Plane
%received - vector of form I+jQ
%ref - reference constellation of form I+jQ
%Note: MPAM/BPSK are one dim. modulations. The same function can be
%applied for these modulations since quadrature is zero (Q=0).
x=[real(received); imag(received)]';%received vec. in cartesian form
y=[real(ref); imag(ref)]';%reference vec. in cartesian form
[idealPoints,indices]= minEuclideanDistance(x,y);
end
function [idealPoints,indices]= minEuclideanDistance(x,y)
%function to compute the pairwise minimum Distance between two
%vectors x and y in p-dimensional signal space and select the
%vectors in y that provides the minimum distances.
% x - a matrix of size mxp
% y - a matrix of size nxp. This acts as a reference against
% which each point in x is compared.
% idealPoints - contain the decoded vector
% indices - indices of the ideal points in reference matrix y
[m,p1] = size(x);[n,p2] = size(y);
if p1~=p2
error('Dimension Mismatch: x and y must have same dimension')
end
X = sum(x.*x,2);
Y = sum(y.*y,2)';
d = X(:,ones(1,n)) + Y(ones(1,m),:) - 2*x*y';%Squared Euclidean Dist.
[~,indices]=min(d,[],2); %Find the minimum value along DIM=2
idealPoints=y(indices,:);
indices=indices.';
end


function [SER] = ser_awgn(EbN0dB,MOD_TYPE,M,COHERENCE)
%Theoretical Symbol Error Rate for various modulations over AWGN
%EbN0dB - list of SNR per bit values
%MOD_TYPE - 'BPSK','PSK','QAM','PAM','FSK'
%M - Modulation level for the chosen modulation
% - For PSK,PAM,FSK M can be any power of 2
% - For QAM M must be even power of 2 (square QAM only)
%Parameter COHERENCE is only applicable for FSK modulation
%COHERENCE = 'coherent' for coherent FSK detection
% = 'noncoherent' for noncoherent FSK detection
gamma_b = 10.^(EbN0dB/10); %SNR per bit in linear scale
gamma_s = log2(M)*gamma_b; %SNR per symbol in linear scale
SER = zeros(size(EbN0dB));
switch lower(MOD_TYPE)
case 'bpsk'
SER=0.5*erfc(sqrt(gamma_b));
case {'psk','mpsk'}
if M==2, %for BPSK
SER=0.5*erfc(sqrt(gamma_b));
else
if M==4, %for QPSK
Q=0.5*erfc(sqrt(gamma_b)); SER=2*Q-Q.^2;
else %for other higher order M-ary PSK
SER=erfc(sqrt(gamma_s)*sin(pi/M));
end
end
case {'qam','mqam'}
SER = 1-(1-(1-1/sqrt(M))*erfc(sqrt(3/2*gamma_s/(M-1)))).^2;
case {'fsk','mfsk'}
if strcmpi(COHERENCE,'coherent'),
for ii=1:length(gamma_s),
fun = @(q) (0.5*erfc((-q - sqrt(2.*gamma_s(ii)))/sqrt(2))).^(M-1).*1/sqrt(2*pi).*exp(-q.^2/2);
SER(ii) = 1-integral(fun,-inf,inf);
end
else %Default compute for noncoherent
for jj=1:length(gamma_s),
summ=0;
for i=1:M-1,
n=M-1; r=i; %for nCr formula
summ=summ+(-1).^(i+1)./(i+1).*prod((n-r+1:n)./(1:r)).*exp(-i./(i+1).*gamma_s(jj));
end
SER(jj)=summ; %Theoretical SER for non-coherent detection
end
end
case {'pam','mpam'}
SER=2*(1-1/M)*0.5*erfc(sqrt(3*gamma_s/(M^2-1)));
otherwise
display 'ser_awgn.m: Invalid modulation (MOD_TYPE) selected'
end
end

function [grayCoded]=dec2gray(decInput)
%convert decimal to Gray code representation
%example: x= [0 1 2 3 4 5 6 7] %decimal
%y=dec2gray(x)
%returns y = [ 0 1 3 2 6 7 5 4] %Gray coded
[rows,cols]=size(decInput);
grayCoded=zeros(rows,cols);
for i=1:rows
for j=1:cols
grayCoded(i,j)=bitxor(bitshift(decInput(i,j),-1),decInput(i,j));
end;end
end
