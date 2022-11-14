clc 
clear all
close all
nacaseries = input('Enter the 4-digit naca series = ');
 c = input('Enter the chord length = ');
 s = num2str(nacaseries);

 % creating points on airfoil
if numel(s)==2
    s1 = str2double(s(1));
 s2 = str2double(s(2));
 m=0;p=0;t=(10*s1+s2)*0.01;
else
    s1 = str2double(s(1));
 s2 = str2double(s(2));
 s3 = str2double(s(3));
 s4 = str2double(s(4));
    m = s1*0.01; p = s2*0.1 ; t = (10*s3+s4)*0.01;
end

n = 250;
for i= n:-1:1
    
    theta = (i-n)*2*pi/n;
    %theta = 2*pi-theta;
   % angle = theta*180/pi
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
xa(n+1)= c ; 
ya(n+1) = 0; 
yc(n+1) = 0;  % trailing edge

% computing control points and panel length
for i = 1:n
    xmid(i) = (xa(i)+xa(i+1))/2;
    ymid(i) = (ya(i)+ya(i+1))/2;
    Sj(i) = sqrt((xa(i+1)-xa(i))^2+(ya(i+1)-ya(i))^2);
end  
%a = 0;
%alpha = pi/180*a;
 % Calculating angles and integrals   
 
 for i = 1:n  
      phi(i) = atan2((ya(i+1)-ya(i)),(xa(i+1)-xa(i)));
      eta(i) = phi(i)+pi/2;
         %phii(i) = phii(i)*180/pi;
        
 end 


     
for i = 1:n 
    for j = 1:n 
         if i~=j
        
        A=-(xmid(i)-xa(j))*cos(phi(j))-(ymid(i)-ya(j))*sin(phi(j));
        B = (xmid(i)-xa(j))^2+(ymid(i)-ya(j))^2;
        C = sin(phi(i)-phi(j));
        D = (ymid(i)-ya(j))*cos(phi(i))-(xmid(i)-xa(j))*sin(phi(i));
        E = sqrt(B-A^2);
        I(i,j) = C/2*log((Sj(j)^2+2*A*Sj(j)+B)/B)+(D-A*C)/E*(atan2((Sj(j)+A),E)-atan2(A,E));
        J(i,j) = (D-A*C)/(2*E)*log((Sj(j)^2+2*A*Sj(j)+B)/B)-C*(atan2((Sj(j)+A),E)-atan2(A,E));
        end
       end   
end

% solving matrix equation
P = zeros(n,n);

for i = 1:n 
    for j = 1:n 
        if i==j
        P(i,j) = pi;
        else 
        P(i,j) = I(i,j);
        end 
    end
end


Q = zeros(n,1); 
U = 1;

for i = 1:n
    Q(i,1) = -2*pi*U*cos(eta(i));
end 

lamda = zeros(n,1);
lamda = mldivide(P,Q);


% calculating tangential velocity accordng to the equation (17)
for i = 1:n
b(i) = 0.0;
end 
for i = 1:n
    for j = 1:n
        if i~=j
        b(i) = b(i)+lamda(j,1)/2/pi*J(i,j);
        end 
    end 
end

  for i = 1:n
        v(i) = U*sin(eta(i))+b(i);
      ut = v(i)/U;
  cp(i) = 1-ut^2; % calculating cp at each control points 
  
  end
 
 
    plot(xa(1:n/2),cp(1:n/2),'-*r')
    set(gca,'Ydir','reverse')
    hold on 
    plot(xa(n/2+1:n),cp(n/2+1:n),'-*b')
    
    hold on 
    plot(xa,ya,'-k')
    xlabel('x')
    ylabel('cp')
    title('cp vs x')
  

    






