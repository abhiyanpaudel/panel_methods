clc 
clear all
close all
nacaseries = input('Enter the 4-digit naca series = ','s');
 c = input('Enter the chord length = ');

    s1 = str2double(nacaseries(1));
 s2 = str2double(nacaseries(2));
 s3 = str2double(nacaseries(3));
 s4 = str2double(nacaseries(4));
    m = s1*0.01; p = s2*0.1 ; t = (10*s3+s4)*0.01;
% end

n = 60;
for i= 1:n
    
    theta = (i-1)*2*pi/n;
    x = 0.5*c*(1+cos(theta));
if(x/c)<p
    yc(i) = m*c/p^2*(2*p*(x/c)-(x/c)^2);
    dydx(i) = (2*m/p^2)* (p-x/c);
    beta(i) = atan(dydx(i));
else
    yc(i) = m*c/(1-p)^2 * ((1-2*p)+2*p*(x/c)-(x/c)^2);
    dydx(i) = (2*m/(1-p)^2)* (p-x/c);
    beta(i) = atan(dydx(i));
end
yt=5*t*c*(0.2969*sqrt(x/c)-0.1260*(x/c)...
    -0.3516*(x/c)^2+0.2843*(x/c)^3-0.1036*(x/c)^4);

% plot(x,yc,'*r')
% hold on

if(i<(0.5*n+1))
    xa(i)=x - yt*sin(beta(i));
    ya(i)=yc(i)+yt*cos(beta(i));
else
    xa(i)=x + yt*sin(beta(i));
    ya(i)=yc(i)-yt*cos(beta(i));
end

end
xa(n+1)= c ; 
ya(n+1) = 0; 
yc(n+1) = 0;  % trailing edge
% plot(xa,ya,'k-')
%hold on
plot(xa,ya,'k -')
 plot(xa(1:0.5*n+1),yc(1:0.5*n+1))
axis equal
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)

 patch(xa,ya,'g')