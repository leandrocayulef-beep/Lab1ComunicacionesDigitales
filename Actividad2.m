fm = 100000;
tm = 1/fm;
ls = 2000;
f_c = 1000;
f_s = 5000;
t_s = 1/f_s;
tau = 0.5*t_s;
d = tau/t_s;
t = (0:ls-1)*tm;
m_t = sin(2*pi*f_c*t);
NFFT = 2^nextpow2(length(m_t));
f = fm*(0:(NFFT/2))/NFFT;
M_f = fft(m_t, NFFT)/length(m_t);
r = floor(t_s/tm);  % Intervalo entre muestras reales
s = floor(tau/tm);  % Duración del pulso en muestras internas
s_nat = zeros(1,length(t));
for i = 1:r:length(m_t)
    s_nat(i:i+s) = 1;   % Pulsos rectangulares
end
s_nat = s_nat(1:length(t));
m_t_nat = m_t .* s_nat;
M_nat_f = fft(m_t_nat, NFFT)/length(m_t_nat);
m_t_inst = zeros(1,length(t));
for i = 1:r:length(m_t)
    m_t_inst(i:i+s) = m_t(i);
end
m_t_inst = m_t_inst(1:length(t));
M_inst_f = fft(m_t_inst, NFFT)/length(m_t_inst);
figure;
plot(f, abs(M_f(1:NFFT/2+1)), 'b', 'LineWidth', 1.2); hold on;
plot(f, abs(M_nat_f(1:NFFT/2+1)), 'r', 'LineWidth', 1.2);
plot(f, abs(M_inst_f(1:NFFT/2+1)), 'g', 'LineWidth', 1.2);
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Magnitud |M(f)|');
legend('Original', 'PAM Natural', 'PAM Instantáneo');
title('Transformada de Fourier: Señal Original, PAM Natural e Instantáneo');

N = 8;  % Número de bits PCM (configurable)
L = 2^N;  % Niveles de cuantización
m_max = max(m_t_inst);
m_min = min(m_t_inst);
delta = (m_max - m_min)/L;
m_q = round((m_t_inst - m_min)/delta);
m_q(m_q > L-1) = L-1;
m_q(m_q < 0) = 0;
m_pcm = m_q * delta + m_min;
error_q = m_t_inst - m_pcm;
figure;
plot(t, m_t, 'b', 'LineWidth', 1.2); hold on;
plot(t, m_t_inst, 'g', 'LineWidth', 1.2);
stairs(t, m_pcm, 'r', 'LineWidth', 1.2);
grid on;
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('Señal original m(t)', 'PAM instantáneo', 'PCM cuantificado');
title(['PCM con N = ', num2str(N), ' bits']);
figure;
stem(t, error_q, 'm', 'filled');
grid on;
xlabel('Tiempo (s)');
ylabel('Error');
title('Error de cuantización (muestras vs PCM)')
