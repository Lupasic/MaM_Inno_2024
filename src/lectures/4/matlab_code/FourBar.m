clc; clear all;
% Parameters
% https://www.youtube.com/watch?v=5U9Dup_ZcGA

A = 2; C = 4; B = 3; D = 3.1;


t = 0:0.05:100;
theta_init = deg2rad(83.9);

ang_speed = 2;
theta = ang_speed*t+theta_init;

P1 = [0;0];
P4 = D*[1;0];

P2 = A*[cos(theta); sin(theta)]; 
E = sqrt(A^2 + D^2 - 2*A*D*cos(theta));
alfa = asin(A*sin(theta)./E);
beta = acos((E.^2 + B^2 - C^2)./(2*E*B));
P3 = [D - B*cos(alfa+beta); B*sin(alfa+beta)];

 gamma = pi - (alfa + beta);

P3_x = P3(1,:);
P3_y = P3(2,:);

P3_vx = diff(P3_x)./diff(t);
P3_vy = diff(P3_y)./diff(t);

P3_v = sqrt(P3_vx.^2 + P3_vy.^2);

for i=1:length(t)
   axis('equal');
   set(gca,'XLim',[-5 8],'YLim',[-3 7]);

%    ani = subplot(2,1,1);
   P1_circle = viscircles(P1',0.05);
   P2_circle = viscircles(P2(:,i)',0.05);
   P3_circle = viscircles(P3(:,i)',0.05);
   P4_circle = viscircles(P4',0.05); 
   
   A_bar = line([P1(1) P2(1,i)],[P1(2) P2(2,i)]);
   B_bar = line([P2(1,i) P3(1,i)],[P2(2,i) P3(2,i)]);
   C_bar = line([P3(1,i) P4(1)],[P3(2,i) P4(2)]);
   

   
   str1 = ['P3: (' num2str(P3(1,i)) ', ' num2str(P3(2,i)) ')'];
   str2 = ['Time elapsed: '  num2str(t(i)) ' s; ' 'element: ' num2str(i) '; a = ' num2str(A) ', C = ' num2str(C) ', B = ' num2str(B) ', D = ' num2str(D)];
   P3_text = text(P3(1,i),P3(2,i)+0.4,str1);
   Time = text(-2,6,str2);
   theta_text = text(P1(1),P1(2)+0.4,['Theta: ' num2str(rad2deg(mod(theta(i),360)))]);
   gamma_text = text(P4(1),P4(2)+0.4,['Gamma: ' num2str(rad2deg(gamma(i)))]);
   pause(1.5);
   if i<length(t)
    delete(P1_circle);
    delete(P2_circle);
%     delete(P3_circle);
    delete(P4_circle);
    delete(A_bar);
    delete(B_bar);
    delete(C_bar);
    delete(P3_text);
    delete(theta_text);
    delete(gamma_text);
    delete(Time);
%     vel = subplot(2,1,2);
%     plot(vel,t(1:i),P3_v(1:i));
%     set(vel,'XLim',[0 10],'YLim',[0 10]);
%     xlabel(vel, 'Time (s)');
%     ylabel(vel, 'Amplitude (m/s');
%     title(vel,'Speed of P3');
%     grid on;
   end

end

