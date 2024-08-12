% MATLAB code to generate plots for symbol error rates

% Define Eb/N0 range
EbN0_dB = 0:2:24;
EbN0 = 10.^(EbN0_dB/10);

% PSK modulation
M_PSK = [2, 4, 8, 16, 32]; % M values for PSK
P_SER_PSK = zeros(length(M_PSK), length(EbN0));
for i = 1:length(M_PSK)
    M = M_PSK(i);
    P_SER_PSK(i, :) = 2*qfunc(sqrt(2*EbN0*sin(pi/M)));
end

% PAM modulation
M_PAM = [2, 4, 8, 16]; % M values for PAM
P_SER_PAM = zeros(length(M_PAM), length(EbN0));
for i = 1:length(M_PAM)
    M = M_PAM(i);
    P_SER_PAM(i, :) = 2*(M-1)/M*qfunc(sqrt(6*log2(M)/(M^2-1)*EbN0));
end

% QAM modulation
M_QAM = [4, 16, 64, 256]; % M values for QAM
P_SER_QAM = zeros(length(M_QAM), length(EbN0));
for i = 1:length(M_QAM)
    M = M_QAM(i);
    P_SER_QAM(i, :) = 4*(1-1/sqrt(M))*qfunc(sqrt(3*log2(M)/(M-1)*EbN0));
end

% MFSK modulation (coherent detection)
M_MFSK = [2, 4, 8, 16, 32]; % M values for MFSK
P_SER_MFSK_coherent = zeros(length(M_MFSK), length(EbN0));
for i = 1:length(M_MFSK)
    M = M_MFSK(i);
    P_SER_MFSK_coherent(i, :) = qfunc(sqrt(2*EbN0*log2(M)/(M-1)));
end

% MFSK modulation (non-coherent detection)
P_SER_MFSK_noncoherent = zeros(length(M_MFSK), length(EbN0));
for i = 1:length(M_MFSK)
    M = M_MFSK(i);
    P_SER_MFSK_noncoherent(i, :) = exp(-EbN0/2).*((1+EbN0).^(M-1)-1)./(M-1);
end

% Plotting
figure;
% PSK plot
subplot(2,3,1);
for i = 1:length(M_PSK)
    semilogy(EbN0_dB, P_SER_PSK(i,:), '-o'); hold on;
end
title('PSK');
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate');
legend(arrayfun(@(x) sprintf('%d-PSK', x), M_PSK, 'UniformOutput', false));
grid on;

% PAM plot
subplot(2,3,2);
for i = 1:length(M_PAM)
    semilogy(EbN0_dB, P_SER_PAM(i,:), '-o'); hold on;
end
title('PAM');
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate');
legend(arrayfun(@(x) sprintf('%d-PAM', x), M_PAM, 'UniformOutput', false));
grid on;

% QAM plot
subplot(2,3,3);
for i = 1:length(M_QAM)
    semilogy(EbN0_dB, P_SER_QAM(i,:), '-o'); hold on;
end
title('QAM');
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate');
legend(arrayfun(@(x) sprintf('%d-QAM', x), M_QAM, 'UniformOutput', false));
grid on;

% MFSK (coherent) plot
subplot(2,3,4);
for i = 1:length(M_MFSK)
    semilogy(EbN0_dB, P_SER_MFSK_coherent(i,:), '-o'); hold on;
end
title('Coherently detected MFSK');
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate');
legend(arrayfun(@(x) sprintf('M=%d', x), M_MFSK, 'UniformOutput', false));
grid on;

% MFSK (non-coherent) plot
subplot(2,3,5);
for i = 1:length(M_MFSK)
    semilogy(EbN0_dB, P_SER_MFSK_noncoherent(i,:), '-o'); hold on;
end
title('Non-coherently detected MFSK');
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate');
legend(arrayfun(@(x) sprintf('M=%d', x), M_MFSK, 'UniformOutput', false));
grid on;

% Adjust subplot positions
set(gcf, 'Position', [100, 100, 1200, 600]);
