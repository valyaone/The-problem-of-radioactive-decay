clear; clc;

U0=10;
lamda=0.5;
N=60;
tau=0.1;
T=[];
U=[];

Un=U0; 
T(1)=0;

   for tn=2:N
   T(tn)=(tn-1)*tau;
   end

%Аналитиическое решение
syms t;
Ut=U0*exp(-lamda*t);

   method=5; %1 - ERK1,2,3,4 // 5 - CROS1
   
if method == 1 %ERK1    
   fu=-lamda*Un;
    for n=1:N
        U(n)=Un;
        Un=Un+tau*fu;
        fu=-lamda*Un;
    end
end

if method == 2 %ERK2
    a2=2/3;
    b1=1/4;
    b2=3/4;
    c1=0;
    c2=0;
    U2(1) = U0;
    for n=2:N
        w1(n-1) = -lamda*U2(n-1);
        
        Uw1=U2+tau*a2*w1;
        
        w2(n-1) = -lamda*Uw1(n-1);
        
        U2(n)=U2(n-1)+tau*(w1(n-1)*b1+b2*w2(n-1));
    end
    U = U2;
end

if method == 3 %ERK3
    a2=1/2;
    a3=3/4; 
    b1=2/9;
    b2=3/9;
    b3=4/9; 
    U2(1) = U0;
    for n=2:N
        w1(n-1) = -lamda*U2(n-1);
        
        Uw1=U2+tau*a2*w1;
        
        w2(n-1) = -lamda*Uw1(n-1);
        
        Uw2=U2+tau*a3*w2;
        
        w3(n-1) = -lamda*Uw2(n-1);
        
        U2(n)=U2(n-1)+tau*(w1(n-1)*b1+b2*w2(n-1)+b3*w3(n-1));
    end
    U = U2;
end

if method == 4 %ERK4
    a2=1/2;
    a3=1/2; 
    a4=1;
    b1=1/6;
    b2=1/3;
    b3=1/3; 
    b4=1/6;
    U2(1) = U0;
    for n=2:N
        w1(n-1) = -lamda*U2(n-1);
        
        Uw1=U2+tau*a2*w1;
        
        w2(n-1) = -lamda*Uw1(n-1);
        
        Uw2=U2+tau*a3*w2;
        
        w3(n-1) = -lamda*Uw2(n-1);
        
        Uw3=U2+tau*a4*w3;
        
        w4(n-1) = -lamda*Uw3(n-1);
        
        U2(n)=U2(n-1)+tau*(w1(n-1)*b1+b2*w2(n-1)+b3*w3(n-1)+b4*w4(n-1));
    end
    U = U2;
end

if method == 5 %CROS1 
    A = (1+1i)/2; %cros1
    %A = 1; %обратная схема эйлера
    %A = 1/2; %схема с полусуммой
    %A = 0; %явная
    
    for n=1:N
    U(n)=Un; 
    w=( - lamda*Un)/( 1-A*tau*(-lamda) );
    Un=Un+tau*real(w);
    end      
end

if method == 6 %CROS2 - надо дописать
    a11 = 0.1 + (sqrt(11)/30)*1i;
    a22 = 0.2 + 0.1i;
    b1 = 0.1941430241155180 - 0.2246898944678803i;
    b2 = 0.8058569758844820 - 0.8870089521907592i;
    c21 = 0.2554708972958462 - 0.2026195833570109i;
    a21 = 0.5617645150714754 - 1.148223341045841i;
          
end

figure(method);
plot(T,U,'r*')
hold on;
ezplot(t,Ut, [0, 10]); 
hold on;