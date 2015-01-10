clc
clear all
close all
%%

T=1;
t11=-3*T:.01:3*T;
t12=t11-T;
t13=t11+T;
a=.9;
p1k=sinc(t11/T).*cos(pi*a*t11/T)./(1-4*a^2*t11.^2/T^2);

tau=.3;
t21=t11+tau;
t22=t11-T+tau;
t23=t11+T+tau;


plot(t11,p1k,'b','LineWidth',2);
hold on
plot(t21,p1k,'k--','LineWidth',2);
legend('Relay 1',' Relay 2')
plot(t12,p1k,'b','LineWidth',2);
plot(t22,p1k,'k--','LineWidth',2);
plot(t23,p1k,'k--','LineWidth',2);
axis([-2 3 -.2 1])
xlabel('Time');
%legend('x1','x2','x3','x3','x5')
set(gca,'XTick',-2:1:2,'YTick',[],'FontSize',16,...
   'FontName','Times New Roman');
box off
grid on
%%
figure
subplot(2,1,1)
plot(t11,p1k,'b','LineWidth',2);
hold on
plot(t12,p1k,'b--','LineWidth',2);
plot(t13,p1k,'b-.','LineWidth',2);
axis([-2 3 -.1 1])
xlabel('Time');
legend('x1','x2','x3')
set(gca,'XTick',-2:1:2,'YTick',[-.1:.3:1],'FontSize',16,...
   'FontName','Times New Roman');
box off
grid on

subplot(2,1,2)
plot(t21,p1k,'k','LineWidth',2);
hold on
plot(t22,p1k,'k--','LineWidth',2);
plot(t23,p1k,'k-.','LineWidth',2);
axis([-2 3 -.1 1])
xlabel('Time');
legend('x4','x5','x6')
set(gca,'XTick',-2:1:2,'YTick',[-.1:.3:1],'FontSize',16,...
   'FontName','Times New Roman');
box off
grid on

