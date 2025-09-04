clc;
close all;
clear all;
fm = 100000;      % Frecuencia de muestreo interna
tm = 1/fm;        % Conversion a 10 microsegundos (mu)
ls = 200;         % Número de muestras
f_c = 1000;       % Frecuencia sinusoidal
f_s = 5000;       % Frecuencia de muestreo real (externa)
t_s = 1/f_s;      % Periodo de muestreo real
tau = 0.5*t_s;    % Duración del pulso de muestreo (tau)
d = tau/t_s;      % Ciclo de trabajo
%% Generando la señal
t = (0:ls-1)*tm;              % Vector de tiempo
m_t = sin(2*pi*f_c*t);        % Señal senoidal de frecuencia f_c
Fourier_original = fft(m_t);  %Fourier
Fourier_original(Fourier_original<0)=0;
%% Cálculos auxiliares para muestreo
r = floor(t_s/tm);  % Intervalo entre muestras reales (en muestras internas)
s = floor(tau/tm);   % Duración del pulso de muestreo (en muestras internas)
%% Muestreo natural
s_nat = zeros(1,length(t)); % Inicia vector de pulsos
for i = 1:length(m_t)
    if mod(i, r) == 0  % Cada r muestras, genera un pulso de duración s
        s_nat(i:i+s) = 1;
    end
end
s_nat = s_nat(1:length(t));   % Asegura que el vector no exceda el tamaño original
m_t_nat = m_t .* s_nat;  % Multiplica la señal por los pulsos → muestreo natural
Fourier_natural = fft(m_t);  %Fourier
Fourier_natural(Fourier_natural<0)=0;
%% Muestreo instantáneo
m_t_inst = zeros(1,length(t));  % Inicia la señal instantanea
for i = 1:length(m_t)
    if mod(i, r) == 0  % Cada r muestras, copia el valor de la senial
        m_t_inst(i:i+s) = m_t(i);
    end
end
m_t_inst = m_t_inst(1:length(t));  % Asegura que el vector no exceda el tamanio original
%% Figura 1: Comparación de muestreos
figure;
plot(t, m_t, 'b');   % Señal original
hold on;
plot(t, m_t_nat, 'r');  % Señal natural
plot(t, m_t_inst, 'g');  % Señal instantanea
grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('Original', 'Natural', 'Instantáneo');
title('Comparación de Muestreos: Natural vs. Instantáneo');
%% Cuantificación PCM de la señal muestreada instantáneamente
bit_depth = 8;                         % Profundidad de bits (8 bits → 256 niveles)
pcm_levels = 2^bit_depth;              % Número total de niveles PCM
pcm_signal_inst = round((m_t_inst + 1) * (pcm_levels - 1) / 2);  
% Convierte la señal de rango [-1,1] a [0,255] y redondea
%% Normalización para visualización
m_t_norm = (m_t - min(m_t)) / (max(m_t) - min(m_t));                     % Señal original normalizada
m_t_inst_norm = (m_t_inst - min(m_t_inst)) / (max(m_t_inst) - min(m_t_inst)); % Señal muestreada normalizada
pcm_signal_inst_norm = (pcm_signal_inst - min(pcm_signal_inst)) / (max(pcm_signal_inst) - min(pcm_signal_inst)); 
% Señal cuantificada normalizada
%% Cálculo del error de cuantización
quantization_error_inst = m_t_inst - ((2 * pcm_signal_inst / (pcm_levels - 1)) - 1);  
