clc;
close all;
clear all;
fm = 100000;      % Frecuencia de muestreo interna
tm = 1/fm;        % Periodo interno
ls = 2000;        % Número de muestras (aumenté para mejor resolución espectral)
f_c = 1000;       % Frecuencia sinusoidal
f_s = 5000;       % Frecuencia de muestreo real (externa)
t_s = 1/f_s;      % Periodo de muestreo real
tau = 0.5*t_s;    % Duración del pulso de muestreo (tau)
d = tau/t_s;      % Ciclo de trabajo
%% Generando la señal
t = (0:ls-1)*tm;              % Vector de tiempo
m_t = sin(2*pi*f_c*t);        % Señal senoidal
%% FFT de la señal original
NFFT = 2^nextpow2(length(m_t));
f = fm*(0:(NFFT/2))/NFFT;     % Vector de frecuencias (positivas)
M_f = fft(m_t, NFFT)/length(m_t);
%% Cálculos auxiliares para muestreo
r = floor(t_s/tm);  % Intervalo entre muestras reales
s = floor(tau/tm);  % Duración del pulso en muestras internas
%% Muestreo natural (PAM natural)
s_nat = zeros(1,length(t));
for i = 1:r:length(m_t)
    s_nat(i:i+s) = 1;   % Pulsos rectangulares
end
s_nat = s_nat(1:length(t));
m_t_nat = m_t .* s_nat;
%% FFT PAM natural
M_nat_f = fft(m_t_nat, NFFT)/length(m_t_nat);
%% Muestreo instantáneo (PAM instantáneo)
m_t_inst = zeros(1,length(t));
for i = 1:r:length(m_t)
    m_t_inst(i:i+s) = m_t(i);  % Mantiene amplitud constante en cada pulso
end
m_t_inst = m_t_inst(1:length(t));
%% FFT PAM instantáneo
M_inst_f = fft(m_t_inst, NFFT)/length(m_t_inst);
%% Figura: espectros de Fourier
figure;
plot(f, abs(M_f(1:NFFT/2+1)), 'b', 'LineWidth', 1.2); hold on;
plot(f, abs(M_nat_f(1:NFFT/2+1)), 'r', 'LineWidth', 1.2);
plot(f, abs(M_inst_f(1:NFFT/2+1)), 'g', 'LineWidth', 1.2);
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Magnitud |M(f)|');
legend('Original', 'PAM Natural', 'PAM Instantáneo');
title('Transformada de Fourier: Señal Original, PAM Natural e Instantáneo');
