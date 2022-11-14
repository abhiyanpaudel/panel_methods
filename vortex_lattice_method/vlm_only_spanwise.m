% This program calculates the cl of a planar swept wing with panels divided in 
% spanwise direction using vortex lattice method
% written by Abhiyan Paudel
% Date : 12/01/2017

clc
clear all
close all 

a = input('Enter the angle of attack = ');
U = input('Enter the freestream velocity = ');
num = input('Enter the number of panels = ');
b = 1;  % wing span 

% Creating Panels
xa_rt = 0;
ya_rt = 0;

xb_rt = 0.2*b;
yb_rt = 0;

xc_rt = 0.7*b;
yc_rt = 0.5*b;

xd_rt = 0.5*b;
yd_rt = 0.5*b;


plot(xa_rt,ya_rt,'ro',xb_rt,yb_rt,'ro',xc_rt,yc_rt,'ro',xd_rt,yd_rt,'ro')
v1 = [xa_rt,xb_rt];
v2 = [ya_rt,yb_rt];
line(v1,v2)

v1 = [xb_rt,xc_rt];
v2 = [yb_rt,yc_rt];
line(v1,v2)

v1 = [xc_rt,xd_rt];
v2 = [yc_rt,yd_rt];
line(v1,v2)

v1 = [xd_rt,xa_rt];
v2 = [yd_rt,ya_rt];
line(v1,v2)

hold on

dx_rt = (xd_rt-xa_rt)/num;
dy_rt = (yd_rt-ya_rt)/num;

% Creating Collocation points 
for i = 1:num
    xtop_rt(i) = (i-1)*dx_rt+0.2*b/4;
    ytop_rt(i) = (i-1)*dy_rt;
    
end
xtop_rt(num+1) = xd_rt+0.2*b/4;
ytop_rt(num+1) = yd_rt;

plot(xtop_rt,ytop_rt,'*r')

hold on 

for i = 1:num
    xbottom_rt(i) = (i-1)*dx_rt+0.2*b;
    ybottom_rt(i) = (i-1)*dy_rt;
end
xbottom_rt(num+1) = xd_rt+0.2*b;
ybottom_rt(num+1) = yd_rt;


  plot(xbottom_rt,ybottom_rt,'r*')
  hold on
  
  for i = 1:num
  x_ctrl_rt(i) = (xtop_rt(i)+xtop_rt(i+1))/2+0.2*b/2;
  y_ctrl_rt(i) = (ytop_rt(i)+ytop_rt(i+1))/2;
  end 
  plot(x_ctrl_rt,y_ctrl_rt,'*b')
  
  hold on

for i = 1:num
    v1 = [xtop_rt(i),xtop_rt(i+1)];
    v2 = [ytop_rt(i),ytop_rt(i+1)];
line(v1,v2)
end 
hold on

for i = 1:num
    v1 = [xtop_rt(i),xbottom_rt(i)];
    v2 = [ytop_rt(i),ybottom_rt(i)];
line(v1,v2)
end  

% Calculate downwash velocity at mth panel induced by horse shoe vortex of 
% n panels of the starboard wing


for m = 1:num
    for n = 1:num
        x_m1n = x_ctrl_rt(m)-xtop_rt(n);
        y_m1n = y_ctrl_rt(m)-ytop_rt(n);
        x_m2n = x_ctrl_rt(m)-xtop_rt(n+1);
        y_m2n = y_ctrl_rt(m)-ytop_rt(n+1);
        x_2n1n = xtop_rt(n+1)-xtop_rt(n);
        y_2n1n = ytop_rt(n+1)-ytop_rt(n);
        d_m1n = sqrt(x_m1n^2+y_m1n^2);
        d_m2n = sqrt(x_m2n^2+y_m2n^2);
        first_first_term = 1/(x_m1n*y_m2n-x_m2n*y_m1n);
        first_sec_term = ((x_2n1n*x_m1n+y_2n1n*y_m1n)/d_m1n)-((x_2n1n*x_m2n+y_2n1n*y_m2n)/d_m2n);
        first_term = first_first_term*first_sec_term;
        sec_term = (-1/y_m1n)*(1+(x_m1n/d_m1n));
        third_term = (-1/y_m2n)*(1+(x_m2n/d_m2n));
        w_s(m,n) = first_term+sec_term-third_term;     % downwash at starboard wing
    end 
end 

% Calculate downwash velocity at mth panelof startboard wing induced by horse 
% shoe vortex of nth panels of the port wing

% Creating Panels on port wing
xa_lt = xa_rt;
ya_lt = -ya_rt;

xb_lt = xb_rt;
yb_lt = -yb_rt;

xc_lt = xc_rt;
yc_lt =-yc_rt;

xd_lt = xd_rt;
yd_lt = -yd_rt;


plot(xa_lt,ya_lt,'ro',xb_lt,yb_lt,'ro',xc_lt,yc_lt,'ro',xd_lt,yd_lt,'ro')
v1 = [xa_lt,xb_lt];
v2 = [ya_lt,yb_lt];
line(v1,v2)

v1 = [xb_lt,xc_lt];
v2 = [yb_lt,yc_lt];
line(v1,v2)

v1 = [xc_lt,xd_lt];
v2 = [yc_lt,yd_lt];
line(v1,v2)

v1 = [xd_lt,xa_lt];
v2 = [yd_lt,ya_lt];
line(v1,v2)

hold on
% 
% dx_lt = (xd_lt-xa_lt)/num;
% dy_lt = (yd_lt-ya_lt)/num;

for i = 1:num
    xtop_lt(i) = xtop_rt(i);
    ytop_lt(i) = -ytop_rt(i);
    
end
xtop_lt(num+1) = xtop_rt(num+1);
ytop_lt(num+1) = -ytop_rt(num+1);

plot(xtop_lt,ytop_lt,'*r')

hold on 

for i = 1:num
    xbottom_lt(i) = xbottom_rt(i);
    ybottom_lt(i) = -ybottom_rt(i);
end
xbottom_lt(num+1) = xbottom_rt(num+1);
ybottom_lt(num+1) = -ybottom_rt(num+1);


  plot(xbottom_lt,ybottom_lt,'r*')
  hold on
  
  for i = 1:num
  x_ctrl_lt(i) = x_ctrl_rt(i);
  y_ctrl_lt(i) = -y_ctrl_rt(i);
  end 
  plot(x_ctrl_lt,y_ctrl_lt,'*b')
  
  hold on

for i = 1:num
    v1 = [xtop_lt(i),xtop_lt(i+1)];
    v2 = [ytop_lt(i),ytop_lt(i+1)];
line(v1,v2)
end 
hold on

for i = 1:num
    v1 = [xtop_lt(i),xbottom_lt(i)];
    v2 = [ytop_lt(i),ybottom_lt(i)];
line(v1,v2)
end  



for m = 1:num
    for n = 1:num
        x_m1n = x_ctrl_rt(m)-xtop_lt(n+1);
        y_m1n = y_ctrl_rt(m)-ytop_lt(n+1);
        x_m2n = x_ctrl_rt(m)-xtop_lt(n);
        y_m2n = y_ctrl_rt(m)-ytop_lt(n);
        x_2n1n = xtop_lt(n)-xtop_lt(n+1);
        y_2n1n = ytop_lt(n)-ytop_lt(n+1);
        d_m1n = sqrt(x_m1n^2+y_m1n^2);
        d_m2n = sqrt(x_m2n^2+y_m2n^2);
        first_first_term = 1/(x_m1n*y_m2n-x_m2n*y_m1n);
        first_sec_term = ((x_2n1n*x_m1n+y_2n1n*y_m1n)/d_m1n)-((x_2n1n*x_m2n+y_2n1n*y_m2n)/d_m2n);
        first_term = first_first_term*first_sec_term;
        sec_term = (-1/y_m1n)*(1+(x_m1n/d_m1n));
        third_term = (-1/y_m2n)*(1+(x_m2n/d_m2n));
        w_p(m,n) = first_term+sec_term-third_term;     % downwash at starboard wing
    end 
end 


for m = 1:num
    for n = 1:num
        w(m,n) = (w_s(m,n)+w_p(m,n))/(4*pi); % net downwash at starboard wing
    end 
end

alpha = a*pi/180;


for i = 1:num
    rhs(i) = -U*sin(alpha);
end 


rho = 1.225;
% calculating circulation by solving linear equations

g = mldivide(w,rhs');


gamma = 0;
for i = 1:num
    gamma = gamma + g(i)*dy_rt;
end 

L = 2*rho*U*gamma;

c = 0.2;

q = 0.5*rho*U^2;
S = b*c;
cl = L/(q*S)

    


