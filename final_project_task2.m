% Final Project - Task 2
% Andrea Senacheribbe s224178 

clear variables
clc
close all

% signals and filter definitions

% params of the filter
wc=20;
dt=0.005;

t_max=1;
N_tot=t_max/dt; %number of samples
t=(0:(N_tot-1))*dt;

ha= wc*exp(-wc*t)+2*wc*(sqrt(3)/3)*exp(-wc*t/2).*cos(-wc*(sqrt(3)/2)*t+5*pi/6);

%defining the input signal
x=zeros(N_tot,1);
x(1)=1; %discrete delta

y=zeros(N_tot,1);
k=zeros(N_tot,1);

%constants of finite difference equation
A=(wc*dt)/(2+wc*dt);
B=(2-wc*dt)/(2+wc*dt);
C=((wc*dt)^2)/(4+2*wc*dt+(wc*dt)^2);
D=(8-2*(wc*dt)^2)/(4+2*wc*dt+(wc*dt)^2);
E=(-4+2*wc*dt-(wc*dt)^2)/(4+2*wc*dt+(wc*dt)^2);

% computing the output of H1
k(1)=A*x(1);
for i=2:N_tot
    k(i)=A*(x(i)+x(i-1))+B*k(i-1);
end

% computing the output of H2
y(1)=C*k(1);
y(2)=C*(k(2)+2*k(1))+D*y(1);
for i=3:N_tot
     y(i)=C*(k(i)+2*k(i-1)+k(i-2))+D*y(i-1)+E*y(i-2);
end

figure('PaperOrientation','landscape')
stem(t, y/dt), hold on, grid on %plotting discrete impulse response
plot(t,ha) % vs continuos
title('Impulse response of the filter'), legend('discrete', 'continuous')
xlabel('t (sec)'), ylabel('h[t/\Deltat], h_a(t)')
print('latex/graphics/task2/impulse_resp','-dpdf')


% frequency response

%defining discrete transfer function
H1 = @(z) (wc*dt*(1+z.^-1))./((z.^-1)*(-2+wc*dt)+2+wc*dt);
H2 = @(z) (wc*dt)^2*(1+z.^-2+2*z.^-1)./((z.^-2)*(4+(wc*dt)^2-2*wc*dt)+(z.^-1)*(-8+2*(wc*dt)^2)+4+(wc*dt)^2+2*wc*dt);

Haf= @(f) (wc^6)./(wc^6+(2*pi*f).^6); %square modulus of frequency response

df=0.01;
f=-(1/(dt)):df:(1/(dt));
dtft=abs(H1(exp(1j*2*pi*f*dt)).*H2(exp(1j*2*pi*f*dt)));

figure('PaperOrientation','landscape')
plot(f, dtft), hold on, grid on
plot(f, sqrt(Haf(f)))
title('Modulus of frequency response'), legend('DTFT', 'analog')
xlabel('f (Hz)'), ylabel('|DTFT[f]|, |H_a(f)|')
axis([-Inf, Inf, -Inf, 1.05])
print('latex/graphics/task2/frequency_resp_lin','-dpdf')


figure('PaperOrientation','landscape')
semilogx(f, 20*log10(dtft)), hold on, grid on
semilogx(f, 10*log10(Haf(f)))
title('Modulus of frequency response (logaritmic scale)'), legend('DTFT', 'analog')
xlabel('log_{10}(f/f0)'), ylabel('20log_{10}|DTFT[f]|, 20log_{10}|H_a(f)|')
axis([-Inf, 1/(2*dt), 10*log10(Haf(1/(2*dt))),3])
print('latex/graphics/task2/frequency_resp_log','-dpdf')

% FFT

f=(-N_tot/2:N_tot/2-1)/(N_tot*dt);
fft_y=abs(fftshift(fft(y)));

figure('PaperOrientation','landscape')
plot(f, fft_y), grid on, hold on
plot(f, sqrt(Haf(f)))
title('Modulus of FFT and analog frequency response'), legend('FFT', 'analog')
xlabel('f (Hz)'), ylabel('|FFT[f]|, |H_a(f)|')
axis([-Inf, Inf, -Inf, 1.05])
print('latex/graphics/task2/fft','-dpdf')

figure('PaperOrientation','landscape')
plot(f, 20*log10(fft_y)), grid on, hold on
plot(f, 10*log10(Haf(f)))
title('Modulus of FFT and analog frequency response (log scale on y)'), legend('FFT', 'analog')
xlabel('f (Hz)'), ylabel('20log_{10}|FFT[f]|, 20log_{10}|Ha(f)|')
axis([-Inf, Inf, 10*log10(Haf(1/(2*dt))),3])
print('latex/graphics/task2/fft_log','-dpdf')

%energy

ha2= @(t) (wc*exp(-wc*t)+2*wc*(sqrt(3)/3)*exp(-wc*t/2).*cos(-wc*(sqrt(3)/2)*t+5*pi/6)).^2;
energy_the=integral(ha2, 0, Inf)
energy_hn=y'*y/dt
energy_hk=fft_y'*fft_y
