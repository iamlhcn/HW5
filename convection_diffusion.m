clear; close all; clc
n = 51;  %grid has n - 2 interior points per dimension (overlapping)
x = linspace(0,1,n);
t = linspace(0,8,n);
dx = 1/n;

T1 = zeros(1,n);


for i = 1:n
    if (x(i)>= 0 && x(i)<= 0.2)
        T1(i) = 1-(10*x(i)-1).^2;
    else
        T1(i) = 0;
    end
end

T1(n) = 0; %TOP
T1(1) = 0;  %BOTTOM
dt = [1.9223e-05 1.9223e-06];

u = 0.08;

T = T1;

alpha = [u*dt(1)/2*dx u*dt(2)/2*dx]; % alpha = (u*dt)/(2*dx);
beta = [1e-3*dt(1)/dx^2 1e-3*dt(2)/dx^2]; % beta = 1e-3*dt/dx^2;
tMAX = 4;
ta = tMAX;

for ii = 1:length(alpha)
    a = dt(ii);
    for m = 0:a:tMAX
        if (m == 0)
            Told = T1; 
        else
            Told = T;
        end
         for i = 2:n-1
             T(i) = Told(i) - alpha(ii)*(Told(i+1)-Told(i-1)) + beta(ii)*(Told(i+1)-2*Told(i)+Told(i-1)); 
          end
    end
    plot(x,T,'LineWidth', 2)
    hold on    
end

%%%% exact %%%%%
 
y = @(x,ta) 1-(10*(x-u*ta)-1).^2;
T_exact = zeros(1,n);
for i = 1:n
    if (x(i)>= u*ta && x(i)<= 0.2 + u*ta)
        T_exact(i) = y(x(i),ta);
    else
        T_exact(i) = 0;
    end
        
end

plot(x,T_exact, 'r', 'LineWidth', 2)

xlabel('X')
ylabel('T')
title('Explicit Euler at t = 4')
legend({'CFL = 2.e-2','CFL = 2.e-3','Exact'},'Location','northeast')
