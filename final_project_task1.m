% Final Project - Task 1
% Andrea Senacheribbe s224178 

%% signals and filter definitions
clear variables
clc
close all

N=64; % period
N_p=8; % number of periods
N_tot=N*N_p; % total number of samples

n=[0:(N_tot-1)]';

% generating the square wave
DC=0.5; % duty cycle
x_sw=repmat([ones(DC*N,1); -ones((1-DC)*N,1)], N_p, 1);
    % one period of the signal is repeated N_p times
    
% generating the sine wave
x_sin=sqrt(2)*sin(2*pi*n/N);

% filter parameters
a=[-0.9,-0.8,-0.4,0,0.4,0.8,0.9]; 
alpha=1-a;

% outputs of the filter, one column for each val of a
y_sw=zeros(N_tot,length(a)); 
y_sin=zeros(N_tot,length(a));

for i=1:length(a) % iterating for different values of a
    y_sw(1, i)=alpha(i)*x_sw(1);
    y_sin(1, i)=alpha(i)*x_sin(1);
	for j=2:N_tot
        % computing the output using finite difference equations
        y_sw(j, i)=alpha(i)*x_sw(j)+a(i)*y_sw(j-1, i);
        y_sin(j, i)=alpha(i)*x_sin(j)+a(i)*y_sin(j-1, i); 
	end
end

%% plot output square wave
close all
for i=1:length(a)
   figure('PaperOrientation','landscape')
   stem(n,x_sw), hold on, grid on % plotting input
   stem(n, y_sw(:,i)) % vs output of the filter
   title(strcat('Input and output signals for the filter with a=', num2str(a(i)))), legend('input x[n]', 'output y[n]')
   xlabel('n'), ylabel('x[n], y[n]')
   axis([-1 N_tot -2.8 2.8]), pbaspect([2.5 1 1])
   print('-fillpage', strcat('latex/graphics/task1/io_sw_',int2str(i)),'-dpdf')
end

%% plot output sin wave
close all
for i=1:length(a)
   figure('PaperOrientation','landscape')
   stem(n,x_sin), hold on, grid on % plotting input
   stem(n, y_sin(:,i)) % vs output of the filter
   title(strcat('Input and output signals for the filter with a=', num2str(a(i)))), legend('input x[n]', 'output y[n]')
   xlabel('n'), ylabel('x[n], y[n]')
   axis([-1 N_tot -1.5 1.5]), pbaspect([2.5 1 1])
   print('-fillpage', strcat('latex/graphics/task1/io_sin_',int2str(i)),'-dpdf')
end

%% DFT for square wave
close all
X_sw=fft(x_sw(1:N));
for i=1:length(a)
   figure('PaperOrientation','landscape')
   stem(0:N-1,abs(X_sw)), hold on, grid on % dft of input
   stem(0:N-1, abs(fft(y_sw(N*(N_p-1)+1:end,i))))
        % dtf of the output (from last period)
   title(strcat('DFT of input and output signals (a=', num2str(a(i)), ')')), legend('DFT of input x[n]', 'DFT of output y[n]')
   xlabel('k'), ylabel('DFT(x[n]), DFT(y[n])')
   axis([0 N-1 0 45]), pbaspect([2.5 1 1])
   print('-fillpage', strcat('latex/graphics/task1/dft_sw_',int2str(i)),'-dpdf')
end

%% DFT for sin wave
close all
X_sin=fft(x_sin(1:N));

for i=1:length(a)
   figure('PaperOrientation','landscape')
   stem(0:N-1,abs(X_sin)), hold on, grid on % dft of input
   stem(0:N-1, abs(fft(y_sin(N*(N_p-1)+1:end,i))))
        % dtf of the output (from last period)
   title(strcat('DFT of input and output signals (a=', num2str(a(i)), ')')), legend('DFT of input x[n]', 'DFT of output y[n]')
   xlabel('k'), ylabel('DFT(x[n]), DFT(y[n])')
   axis([0 N-1 0 50]), pbaspect([2.5 1 1])
   print('-fillpage', strcat('latex/graphics/task1/dft_sin_',int2str(i)),'-dpdf')
end

%% theoretical sinusoid
close all
H_resp=((1-[-0.8; 0.8])*exp(1j*2*pi/N))./(exp(1j*2*pi/N)-[-0.8; 0.8]);
    % evaluating H(z) (for two val of a) at exp(j 2 pi / N)

% for a = -0.8
figure('PaperOrientation','landscape')
stem(n,y_sin(:,2)), hold on, grid on % plotting filter output
stem(n, abs(H_resp(1))*sqrt(2)*sin(2*pi*n/N+angle(H_resp(1))))
    % vs theoretical output from theoretical sinusoid
title('Comparison of theoretical and simulated output of the filter (a=-0.8)'), legend('simulated output y[n]', 'theoretical output y_{the}[n]')
xlabel('n'), ylabel('y[n], y_{the}[n]')
axis([-1 N*2 -1.5 1.5]), pbaspect([2.5 1 1])
print('-fillpage', 'latex/graphics/task1/theor_sin_2','-dpdf')

% for a = 0.8
figure('PaperOrientation','landscape')
stem(n,y_sin(:,6)), hold on, grid on % plotting filter output 
stem(n, abs(H_resp(2))*sqrt(2)*sin(2*pi*n/N+angle(H_resp(2))))
    % vs theoretical output from theoretical sinusoid
title('Comparison of theoretical and simulated output of the filter (a=0.8)'), legend('simulated output y[n]', 'theoretical output y_{the}[n]')
xlabel('n'), ylabel('y[n], y_{the}[n]')
axis([-1 N*2 -1.5 1.5]), pbaspect([2.5 1 1])
print('-fillpage', 'latex/graphics/task1/theor_sin_6','-dpdf')

%% power
close all
avg_power_sw=zeros(1,length(a));
avg_power_sin=zeros(1,length(a));
for i=1:length(a) % iterating for different values of a
    % evaluating avg power on last period of output signals 
    avg_power_sw(i)=y_sw(N*(N_p-1)+1:end,i)'*y_sw(N*(N_p-1)+1:end,i)/N;
    avg_power_sin(i)=y_sin(N*(N_p-1)+1:end,i)'*y_sin(N*(N_p-1)+1:end,i)/N;
end

figure('PaperOrientation','landscape')
stem(a, avg_power_sw), hold on, grid on
stem(a, avg_power_sin)
title('Comparison of average power for different filter outputs'), legend('square wave', 'sinusoid')
xticks(a), xlabel('a'), ylabel('avg power of y[n]')
pbaspect([2.5 1 1])
print('-fillpage', 'latex/graphics/task1/power','-dpdf')
