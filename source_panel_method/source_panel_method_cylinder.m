clc 
clear all
close all

 % creating points on airfoil


n = 250;
r = 1; % radius of cylinder
for i= n:-1:1
    if i>=n/2
    theta = (i-3*n/2)*2*pi/n;
    else 
        theta = (i-n/2)*2*pi/n;
    end
   
    x(n+1-i) = r*cos(theta);
    y(n+1-i) = r*sin(theta);
end 

x(n+1)= -r ; 
y(n+1) = 0; 

a = 1;
xa(1)=x(n);
ya(1)=y(n);

for i = 1:n
    if mod(i,2)==0
        a = a+1;
        xa(a)=x(i);
        ya(a)=y(i);
    end 
end

% computing control points and panel length
for i = 1:n/2
    xmid(i) = (xa(i)+xa(i+1))/2;
    ymid(i) = (ya(i)+ya(i+1))/2;
    Sj(i) = sqrt((xa(i+1)-xa(i))^2+(ya(i+1)-ya(i))^2);
end  
%a = 0;
%alpha = pi/180*a;
 % Calculating angles and integrals   
 
 for i = 1:n/2  
      phi(i) = atan2((ya(i+1)-ya(i)),(xa(i+1)-xa(i)));
      eta(i) = phi(i)+pi/2;
         %phii(i) = phii(i)*180/pi;
        
 end 


     
for i = 1:n/2 
    for j = 1:n /2
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
P = zeros(n/2,n/2);

for i = 1:n /2
    for j = 1:n/2 
        if i==j
        P(i,j) = pi;
        else 
        P(i,j) = I(i,j);
        end 
    end
end

P
Q = zeros(n/2,1); 
U = 1;

for i = 1:n/2
    Q(i,1) = -2*pi*U*cos(eta(i));
end 

lamda = zeros(n/2,1);
lamda = mldivide(P,Q);


% calculating tangential velocity accordng to the equation (17)
for i = 1:n/2
b(i) = 0.0;
end 
for i = 1:n/2
    for j = 1:n/2
        if i~=j
        b(i) = b(i)+lamda(j,1)/2/pi*J(i,j);
        end 
    end 
end

  for i = 1:n/2
        v(i) = U*sin(eta(i))+b(i);
      ut = v(i)/U;
  cp(i) = 1-ut^2; % calculating cp at each control points 
  
  end
 
  for i = 1:n/2
      theta(i) = (i-1)*2*pi/n/2;
  end
  
plot(theta(1:n/2), cp(1:n/2),'LineWidth', 3);
set(gca,'FontName','Symbol')
set(gca,'XTickLabel','0|  |p/2|  |p|  |3p/2|  |2p')
set(gca,'YTick',-3:1:1)








