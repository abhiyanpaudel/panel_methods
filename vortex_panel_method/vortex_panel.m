% written by Abhiyan Paudel on 07-01-2017
clc 
clear all
close all
nacaseries = input('Enter the 4-digit naca series = ','s');
 c = input('Enter the chord length = ');
 a = input('Enter the angle of attack in degree = ');

 n = 200;
 % creating nodes on airfoil
  s1 = str2num(nacaseries(1));
  s2 = str2num(nacaseries(2));
  s3 = str2num(nacaseries(3:4));
   m = s1*0.01; p = s2*0.1 ; t = s3*0.01;



for i= n:-1:1
    
    theta = (i-n)*2*pi/n;
    x = 0.5*c*(1+cos(theta));
    
if(x/c)<p
    yc(n+1-i) = m*c/p^2*(2*p*(x/c)-(x/c)^2);
    dydx(n+1-i) = (2*m/p^2)* (p-x/c);
    beta(n+1-i) = atan(dydx(n+1-i));
else
    yc(n+1-i) = m*c/(1-p)^2 * ((1-2*p)+2*p*(x/c)-(x/c)^2);
    dydx(n+1-i) = (2*m/(1-p)^2)* (p-x/c);
    beta(n+1-i) = atan(dydx(n+1-i));
end
yt=5*t*c*(0.2969*sqrt(x/c)-0.1260*(x/c)...
    -0.3516*(x/c)^2+0.2843*(x/c)^3-0.1036*(x/c)^4);

if(i<(0.5*n+1))
    xa(n+1-i)=x - yt*sin(beta(n+1-i));
    ya(n+1-i)=yc(n+1-i)+yt*cos(beta(n+1-i));
else
    xa(n+1-i)=x + yt*sin(beta(n+1-i));
    ya(n+1-i)=yc(n+1-i)-yt*cos(beta(n+1-i));
end

end
xa(n+1)= c; 
ya(n+1) = 0; 
yc(n+1) = 0;  % trailing edge


 alpha = a*pi/180;
%This loop calculates the location of the control points 

for i = 1:n
    xmid(i) = (xa(i)+xa(i+1))/2;
    ymid(i) = (ya(i)+ya(i+1))/2;
    Sj(i) = sqrt((xa(i+1)-xa(i))^2+(ya(i+1)-ya(i))^2);%array of panel lengths
    phi(i) = atan2((ya(i+1)-ya(i)),(xa(i+1)-xa(i)));
    rhs(i)= sin(phi(i)-alpha);% RHS of equation 4
end 

% This loop calculates the coefficients to find the matrix Aij for equation
% 6.Variables are as given in report.
for i = 1:n
    for j = 1:n
        if i==j
            cn1(i,j) = -1;
            cn2(i,j) = 1;
            ct1(i,j) = pi/2;
            ct2(i,j) = pi/2;
        else
        A=-(xmid(i)-xa(j))*cos(phi(j))-(ymid(i)-ya(j))*sin(phi(j));
        B = (xmid(i)-xa(j))^2+(ymid(i)-ya(j))^2;
        C = sin(phi(i)-phi(j));
        D = cos(phi(i)-phi(j));
        E = (xmid(i)-xa(j))*sin(phi(j))-(ymid(i)-ya(j))*cos(phi(j));
        F = log(1+(Sj(j)^2+2*A*Sj(j))/B);
        G = atan2((E*Sj(j)),(B+A*Sj(j)));
        P = (xmid(i)-xa(j))*sin(phi(i)-2*phi(j))+(ymid(i)-ya(j))*cos(phi(i)-2*phi(j));
        Q = (xmid(i)-xa(j))*cos(phi(i)-2*phi(j))-(ymid(i)-ya(j))*sin(phi(i)-2*phi(j));
       cn2(i,j)=D+0.5*Q*F/Sj(j)-(A*C+D*E)*G/Sj(j);
       cn1(i,j)=0.5*D*F+C*G-cn2(i,j);
       ct2(i,j)=C+0.5*P*F/Sj(j)+(A*D-C*E)*G/Sj(j);
       ct1(i,j)=0.5*C*F-D*G-ct2(i,j);
        end
    end
end

%This loop calculates the coefficients Anij and Atij used in eq. 6 and
%eq. 8 respectively

for i=1:n
        an(i,1)=cn1(i,1);
        an(i,n+1)=cn2(i,n);
        at(i,1)=ct1(i,1);
        at(i,n+1)=ct2(i,n);
        for j=2:n
            an(i,j)=cn1(i,j)+cn2(i,(j-1));
            at(i,j)=ct1(i,j)+ct2(i,(j-1));
        end
end
    
%kutta condition 

  an(n+1,1)=1;
    an(n+1,n+1)=1;
    rhs(n+1)=0;
    for j=2:n
        an(n+1,j)=0;
    end
        
  
    % calculating circulation density by solving linear equations given by
    % eq. 6
    g = an\rhs';
    
  
    
    % calculating tangential velocity and coefficent of pressure
      for i=1:n
        sum=0;
        for j=1:n+1;
            sum=sum+at(i,j)*g(j);
         end
        v(i) = (cos(phi(i)-alpha)+sum);
        cp(i) = 1-v(i)*v(i);
      end 
      % compute coefficients of lift and drag
      
      cy = 0.0;
      cx = 0.0;
      for i = 1:n
          dx = xa(i+1)-xa(i);
          dy = ya(i+1)-ya(i);
          cy = cy-cp(i)*dx;
          cx = cx+cp(i)*dy;
      end
    
      %  cl and cd 
       cl = cy*cos(alpha)-cx*sin(alpha)
       cd = cy*sin(alpha)+cx*cos(alpha)
      
%    plot(xa,ya,'-k')
%       hold on
    plot(xmid(1:n/2)/c,cp(1:n/2),'-*r')
    set(gca,'Ydir','reverse')
    hold on 
    plot(xmid(n/2+1:n)/c,cp(n/2+1:n),'-*b')
    hold on

  
 
 
    xlabel('x/c')
    ylabel('cp')
    title('cp vs x/c')
    









