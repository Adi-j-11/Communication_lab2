%Eb/N0 Vs SER for M-QAM modulation with receiver impairments
clear all;clc;
%---------Input Fields------------------------
N=100000; %Number of input symbols
EbN0dB = -4:5:24; %Define EbN0dB range for simulation
M=64; %M-QAM modulation order
g=0.9; phi=8; dc_i=1.9; dc_q=1.7;%receiver impairments
%----------------------------------------------
k=log2(M); %Bits per symbol
EsN0dB = 10*log10(k)+EbN0dB; %Converting Eb/N0 to Es/N0
SER1 = zeros(length(EsN0dB),1);%Symbol Error rates (No compensation)
SER2 = SER1;%Symbol Error rates (DC compensation only)
SER3 = SER1;%Symbol Error rates (DC comp & Blind IQ compensation)
SER4 = SER1;%Symbol Error rates (DC comp & Pilot IQ compensation)
d=ceil(M.*rand(1,N));%random data symbols drawn from [1,2,..,M]
[s,ref]=mqam_modulator(M,d);%MQAM symbols & reference constellation
%See section 3.4 on chapter 'Digital Modulators and Demodulators
%- Complex Baseband Equivalent Models' for function definition
for i=1:length(EsN0dB),
r = add_awgn_noise(s,EsN0dB(i)); %see section 4.1 on chapter 4
z=receiver_impairments(r,g,phi,dc_i,dc_q);%add impairments
v=dc_compensation(z); %DC compensation
y3=blind_iq_compensation(v); %blind IQ compensation
[Kest,Pest]=pilot_iq_imb_est(g,phi,dc_i,dc_q);%Pilot based estimation
y4=iqImb_compensation(v,Kest,Pest);%IQ comp. using estimated values
%Enable this section - if you want to plot constellation diagram -
%figure(1);plot(real(z),imag(z),'rO'); hold on;
%plot(real(y4),imag(y4),'b*'); hold off;
%title(['Eb/N0=',num2str(EbN0dB(j)),' (dB)']);pause;
%-------IQ Detectors - defined in section 3.5.4 chapter 3--------
[estTxSymbols_1,dcap_1]= iqOptDetector(z,ref);%No compensation
[estTxSymbols_2,dcap_2]= iqOptDetector(v,ref);%DC compensation only
[estTxSymbols_3,dcap_3]= iqOptDetector(y3,ref);%DC & blind IQ comp.
[estTxSymbols_4,dcap_4]= iqOptDetector(y4,ref);%DC & pilot IQ comp.
%------ Symbol Error Rate Computation-------
SER1(i)=sum((d~=dcap_1))/N; SER2(i)=sum((d~=dcap_2))/N;
SER3(i)=sum((d~=dcap_3))/N; SER4(i)=sum((d~=dcap_4))/N;
end
theoreticalSER = ser_awgn(EbN0dB,'MQAM',M); %theoretical SER
figure(2);
semilogy(EbN0dB,SER1,'r*-'); hold on; % Plot results
semilogy(EbN0dB,SER2,'bO-');
semilogy(EbN0dB,SER3,'g^-'); % Corrected marker to '^'
semilogy(EbN0dB,SER4,'m*-');
semilogy(EbN0dB,theoreticalSER,'k');
legend('No compensation','DC comp only',...
'Sim- DC & blind iq comp','Sim- DC & pilot iq comp','Theoretical');
xlabel('E_b/N_0 (dB)');
ylabel('Symbol Error Rate (P_s)');
title('Probability of Symbol Error 64-QAM signals');


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
end
end
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

function z=receiver_impairments(r,g,phi,dc_i,dc_q)
%Function to add receiver impairments to the IQ branches
%[z]=iq_imbalance(r,g,phi) introduces DC and IQ imbalances between inphase
% and quadrature components of the complex baseband signal r. The model
% parameter g represents gain mismatch between the IQ branches of the receiver
% and parameter 'phi' represents phase error of local oscillator (in degrees).
% DC biases associated with each I,Q path are represented by dc_i and dc_q.
k = iq_imbalance(r,g,phi); %Add IQ imbalance
z = dc_impairment(k,dc_i,dc_q); %Add DC impairment
end

function [z]= iq_imbalance(r,g,phi)
%Function to create IQ imbalance impairment in a complex baseband
% [z]=iq_imbalance(r,g,phi) introduces IQ imbalance and phase error
% signal between the inphase and quadrature components of the
% complex baseband signal r. The model parameter g represents the
% gain mismatch between the IQ branches of the receiver and 'phi'
% represents the phase error of the local oscillator (in degrees).
Ri=real(r); Rq=imag(r);
Zi= Ri; %I branch
Zq= g*(-sin(phi/180*pi)*Ri + cos(phi/180*pi)*Rq);%Q branch crosstalk
z=Zi+1i*Zq;
end

function [y]=dc_impairment(x,dc_i,dc_q)
%Function to create DC impairments in a complex baseband model
% [y]=iq_imbalance(x,dc_i,dc_q) introduces DC imbalance
% between the inphase and quadrature components of the complex
% baseband signal x. The DC biases associated with each I,Q path
% are represented by the paramters dc_i and dc_q
y = x + (dc_i+1i*dc_q);
end

function [v]=dc_compensation(z)
%Function to estimate and remove DC impairments in the IQ branch
% v=dc_compensation(z) removes the estimated DC impairment
iDCest=mean(real(z));%estimated DC on I branch
qDCest=mean(imag(z));%estimated DC on I branch
v=z-(iDCest+1i*qDCest);%remove estimated DCs
end

function y=blind_iq_compensation(z)
%Function to estimate and compensate IQ impairments for the single-
%branch IQ impairment model
% y=blind_iq_compensation(z) estimates and compensates IQ imbalance
I=real(z);
Q=imag(z);
theta1=(-1)*mean(sign(I).*Q);
theta2=mean(abs(I));
theta3=mean(abs(Q));
c1=theta1/theta2;
c2=sqrt((theta3^2-theta1^2)/theta2^2);
yI = I;
yQ = (c1*I+Q)/c2;
y= (yI +1i*yQ);
end

function [Kest,Pest]=pilot_iq_imb_est(g,phi,dc_i,dc_q)
%Length 64 - Long Preamble as defined in the IEEE 802.11a
preamble_freqDomain = [0,0,0,0,0,0,1,1,-1,-1,1,1,-1,1,...
-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,...
0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,...
-1,-1,1,-1,1,-1,1,1,1,1,0,0,0,0,0];%freq. domain representation
preamble=ifft(preamble_freqDomain,64);%time domain representation
%send known preamble through DC & IQ imbalance model and estimate it
r=receiver_impairments(preamble,g,phi,dc_i,dc_q);
z=dc_compensation(r); %remove DC imb. before IQ imbalance estimation
%IQ imbalance estimation
I=real(z); Q=imag(z);
Kest = sqrt(sum((Q.*Q))./sum(I.*I)); %estimate gain imbalance
Pest = sum(I.*Q)./sum(I.*I); %estimate phase mismatch
end

function y=iqImb_compensation(d,Kest,Pest)
%Function to compensate IQ imbalance during the data transmission
% y=iqImb_compensation(d,Kest,Pest) compensates the IQ imbalance
% present at the received complex signal d at the baseband
% processor. The IQ compensation is performed using the gain
% imbalance (Kest) and phase error (Pest) parameters that are
% estimated during the preamble transmission.
I=real(d); Q=imag(d);
wi= I;
wq = (Q - Pest*I)/sqrt(1-Pest^2)/Kest;
y = wi + 1i*wq;
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