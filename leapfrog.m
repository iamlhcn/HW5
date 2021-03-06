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

T1(n) = 0;
T1(1) = 0;  

dt = [1e-2 1e-3];

u = 0.08;

alpha = [0.0204 0.0020]; % alpha = (u*dt)/(2*dx);

tMAX = 8; % maximum Time 0, 4, 8
ta = tMAX;

T_1 = zeros(1,n);
T = T1;
Told = T1;


for ii = 1:length(alpha)
    a = dt(ii);
    for i = 2:n-1
        T_1(i) = T1(i) - alpha(ii)*(T1(i+1)-T1(i-1)); 
    end

    for m = 0:a:tMAX
        if (mod(m,dt) == 0)
              for i = 2:n-1
                  T(i) = Told(i) - 2*alpha(ii)*(Told(i+1)-Told(i-1)); 
              end
              Told = T;
        else
            for i = 2:n-1
                  T(i) = T_1(i) - 2*alpha(ii)*(T_1(i+1)-T_1(i-1)); 
            end
             T_1 = T;
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
title('Leapfrog at t = 8')
legend({'CFL = 2.e-2','CFL = 2.e-3','Exact'},'Location','northwest')

