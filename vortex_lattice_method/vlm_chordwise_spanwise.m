%   This program calculates Cl of planar swept wing with panels divided in
%   chordwise and spanwise directions using vortex lattice method
%   written by Abhiyan Paudel
%   Date : 11/01/2017

clc
clear all
close all

a = input('Enter the angle of attack = ');
U = input('Enter the freestream velocity = ');
N = input('Enter the number of panels in spanwise direction = ');
M = input('Enter the number of panels in chordwise direction = ');


b = 1; % wing span
LX=0.2*b;
LY=0.5*b;

dx=LX/M;
dy=LY/N;

% creating panels on right wing
for i=1:M+1
    for j=1:N+1
        x0_rt(i,j)=(i-1)*dx;
        y_rt(i,j)=(j-1)*dy;
        x_rt(i,j)= x0_rt(i,j)+y_rt(i,j);

    end
end

% computing collocation points
for i = 1:M
    for j = 1:N+1
        xtop_rt(i,j) = x_rt(i,j)+LX/4/M;
        ytop_rt(i,j) =  y_rt(i,j);
    end
end

%computing control points on each panel
k = 0;
for i = 1:M
    for j = 1:N
        k = k+1;
        xctrl_rt(k) = (xtop_rt(i,j)+xtop_rt(i,j+1))/2+LX*0.5/M;
        yctrl_rt(k) = (y_rt(i,j)+y_rt(i,j+1))*0.5;
    end
end

%creating panels on left wing
for i = 1:M+1
    for j = 1:N+1
        x_lt(i,j) = x_rt(i,j);
        y_lt(i,j) = -y_rt(i,j);

    end
end

%computing collocation points
for i = 1:M
    for j = 1:N+1
        xtop_lt(i,j) = xtop_rt(i,j);
        ytop_lt(i,j) =  -ytop_rt(i,j);
    end
end

%computing control points
for i = 1:k

    xctrl_lt(i) = xctrl_rt(i);
    yctrl_lt(i) = -yctrl_rt(i);
end


% calculate downwash at mth panel of starboard wing induced by horse-shoe
% vortex of n panels of the starboard wing
for l = 1:k
    m = 0;
    for i = 1:M
        for j = 1:N
            m = m+1;
            x_m1n = xctrl_rt(l)-xtop_rt(i,j);
            y_m1n = yctrl_rt(l)-ytop_rt(i,j);
            x_m2n = xctrl_rt(l)-xtop_rt(i,j+1);
            y_m2n = yctrl_rt(l)-ytop_rt(i,j+1);
            x_2n1n = xtop_rt(i,j+1)-xtop_rt(i,j);
            y_2n1n = ytop_rt(i,j+1)-ytop_rt(i,j);
            d_m1n = sqrt(x_m1n^2+y_m1n^2);
            d_m2n = sqrt(x_m2n^2+y_m2n^2);
            first_first_term = 1/(x_m1n*y_m2n-x_m2n*y_m1n);
            first_sec_term = ((x_2n1n*x_m1n+y_2n1n*y_m1n)/d_m1n)-((x_2n1n*x_m2n+y_2n1n*y_m2n)/d_m2n);
            first_term = first_first_term*first_sec_term;
            sec_term = (-1/y_m1n)*(1+(x_m1n/d_m1n));
            third_term = (-1/y_m2n)*(1+(x_m2n/d_m2n));
            w_s(l,m) = first_term+sec_term-third_term;     % downwash at starboard wing
        end
    end
end


% calculate downwash at mth panel of starboard wing induced by horse-shoe
% vortex of n panels of the port wing
for l = 1:k
    m = 0;
    for i = 1:M
        for j = 1:N
            m = m+1;
            x_m1n = xctrl_rt(l)-xtop_lt(i,j+1);
            y_m1n = yctrl_rt(l)-ytop_lt(i,j+1);
            x_m2n = xctrl_rt(l)-xtop_lt(i,j);
            y_m2n = yctrl_rt(l)-ytop_lt(i,j);
            x_2n1n = xtop_lt(i,j)-xtop_lt(i,j+1);
            y_2n1n = ytop_lt(i,j)-ytop_lt(i,j+1);
            d_m1n = sqrt(x_m1n^2+y_m1n^2);
            d_m2n = sqrt(x_m2n^2+y_m2n^2);
            first_first_term = 1/(x_m1n*y_m2n-x_m2n*y_m1n);
            first_sec_term = ((x_2n1n*x_m1n+y_2n1n*y_m1n)/d_m1n)-((x_2n1n*x_m2n+y_2n1n*y_m2n)/d_m2n);
            first_term = first_first_term*first_sec_term;
            sec_term = (-1/y_m1n)*(1+(x_m1n/d_m1n));
            third_term = (-1/y_m2n)*(1+(x_m2n/d_m2n));
            w_p(l,m) = first_term+sec_term-third_term;     % downwash at starboard wing
        end
    end
end


for i = 1:k
    for  j= 1:m
        w(i,j) = (w_s(i,j)+w_p(i,j))/(4*pi); % net downwash at starboard wing
    end
end

alpha = a*pi/180;


for i = 1:k
    rhs(i) = -U*sin(alpha);
end


rho = 1.225;
% calculating circulation by solving linear equations

g = mldivide(w,rhs');


gamma = 0;
for i = 1:k
    gamma = gamma + g(i)*dy;
end

L = 2*rho*U*gamma;

c = 0.2;

q = 0.5*rho*U^2;
S = b*c;
cl = L/(q*S)

plot(x_rt,y_rt,'*-r',x_rt',y_rt','-*r')

hold on

plot(xtop_rt',ytop_rt','*-b')

hold on
plot(xctrl_rt',yctrl_rt','*g')
hold on
plot(x_lt,y_lt,'*-r',x_lt',y_lt','-*r')

hold on

plot(xtop_lt',ytop_lt','*-b')

hold on
plot(xctrl_lt',yctrl_lt','*g')

grid on
xlabel('x')
ylabel('y')


