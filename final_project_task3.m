% Final Project - Task 3
% Andrea Senacheribbe s224178 

clear variables
clc
close all

% params of the filter
wc=20;
dt=0.005;

N_tot=10000; %number of samples

N0=[1 1e-5 1e-10 1e-15 1e-20];
N_real=5; %number of realisation for each N0

y=zeros(N_tot,1);
k=zeros(N_tot,1);
variances=zeros(N_real,length(N0)); %variances matrix, each column corresponds to a value of N0, each row is a different realisation

%constants of finite difference equation
A=(wc*dt)/(2+wc*dt);
B=(2-wc*dt)/(2+wc*dt);
C=((wc*dt)^2)/(4+2*wc*dt+(wc*dt)^2);
D=(8-2*(wc*dt)^2)/(4+2*wc*dt+(wc*dt)^2);
E=(-4+2*wc*dt-(wc*dt)^2)/(4+2*wc*dt+(wc*dt)^2);

for i_n0=1:length(N0) %iterating for all values of N0
    for i_real=1:N_real %N_real times for each N0
        
        x=randn(N_tot, 1)*sqrt(N0(i_n0)/(2*dt)); %definition of the WGN
        
        %Butterworth filter
        k(1)=A*x(1);
        for i=2:N_tot % computing the output of H1
            k(i)=A*(x(i)+x(i-1))+B*k(i-1);
        end

        y(1)=C*k(1);
        y(2)=C*(k(2)+2*k(1))+D*y(1);
        for i=3:N_tot % computing the output of H2
             y(i)=C*(k(i)+2*k(i-1)+k(i-2))+D*y(i-1)+E*y(i-2);
        end
        
        variances(i_real, i_n0)=var(y); %computing variances
        
        if (i_n0==1)&&(i_real==1)
            %plotting the histogram for one realisation, N0=1
            figure('PaperOrientation','landscape')
            histogram(y, 'Normalization','probability')
            title('Histogram of distribution of the output signal')
            xlabel('y[n]'), ylabel('probability')
            axis([-max(abs(y)), max(abs(y)), -Inf, Inf]) 
            print('-fillpage','latex/graphics/task3/histogram','-dpdf')
        end

    end
end

figure('PaperOrientation','landscape')
loglog(kron(N0,ones(1,N_real)), variances(:), '*'), hold on %plotting the measured variances
loglog(N0, N0*10/3, 'or') %vs the theoretical ones
title('Variances of the output of the filter'), legend('measured', 'theoretical')
xlabel('N_0'), ylabel('\sigma^2(y[n])')
print('latex/graphics/task3/variances','-dpdf')
axis([0.9, 1.1, 2.8 , 3.8]) 
print('latex/graphics/task3/variances_zoom','-dpdf')


N_tot=1000; %reducing the number of samples
n=0:(N_tot-1);

y=zeros(N_tot,1);
k=zeros(N_tot,1);

for i_real=1:2
        x=randn(N_tot, 1)*sqrt(1/(2*dt)); %definition of the WGN
        
        %Butterworth filter
        k(1)=A*x(1);
        for i=2:N_tot % computing the output of H1
            k(i)=A*(x(i)+x(i-1))+B*k(i-1);
        end

        y(1)=C*k(1);
        y(2)=C*(k(2)+2*k(1))+D*y(1);
        for i=3:N_tot % computing the output of H2
             y(i)=C*(k(i)+2*k(i-1)+k(i-2))+D*y(i-1)+E*y(i-2);
        end

        figure('PaperOrientation','landscape')
        plot(n, x), hold on, grid on %plotting the input
        plot(n, y) % vs output
        title(strcat('Input and output of the filter (realisation #', int2str(i_real), ')')), legend('input', 'output')
        xlabel('n'), ylabel('x[n], y[n]')
        axis([-Inf, Inf, -35, 35]), pbaspect([2.5 1 1])
        print('-fillpage', strcat('latex/graphics/task3/io_',int2str(i_real)),'-dpdf')
end
