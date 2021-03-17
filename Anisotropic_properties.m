clear;
close;
clc;

load fre_90_v.mat Ang v;
load fre_90_inQ.mat Ang inQ;
load fra_intersect_pi_2 ang v1;
v1t=v;
inQ1=inQ;

load fre_60_v.mat Ang v;
load fre_60_inQ.mat Ang inQ;
load fra_intersect_pi_3 ang v2;
v2t=v;
inQ2=inQ;

load fre_30_v.mat Ang v;
load fre_30_inQ.mat Ang inQ;
load fra_intersect_pi_6 ang v3;
v3t=v;
inQ3=inQ;

figure(1)
plot(Ang,v1t(:,1)','r','linewidth',1.3);
hold on;
plot(Ang,v1t(:,2)','g','linewidth',1.3);
hold on;
plot(Ang,v1t(:,3)','b','linewidth',1.3);
hold on;
plot(ang,v1(:,1)','r--','linewidth',1.3);
hold on;
plot(ang,v1(:,2)','g--','linewidth',1.3);
hold on;
plot(ang,v1(:,3)','b--','linewidth',1.3);

xlabel('Incidence angle','FontSize',12);
ylabel('{\itv_{p}}(m/s)','FontSize',12);
legend('{\itf} = 0.001 Hz','{\itf} = 10 Hz','{\itf} = 800 Hz');
xlim([0,360]);
set(gca,'XTick',[0 60 120 180 240 300 360]);


figure(2)
plot(Ang,inQ1(:,1),'r','linewidth',1.3);
hold on;
plot(Ang,inQ1(:,2),'g','linewidth',1.3);
hold on;
plot(Ang,inQ1(:,3),'b','linewidth',1.3);

xlabel('Incidence angle','FontSize',12);
ylabel('1/{\itQ_{p}}','FontSize',12);
legend('{\itf} = 0.1 Hz','{\itf} = 80 Hz','{\itf} = 2500 Hz');
xlim([0,360]);
set(gca,'XTick',[0 60 120 180 240 300 360]);

figure(3)
plot(Ang,v2t(:,1)','r','linewidth',1.3);
hold on;
plot(Ang,v2t(:,2)','g','linewidth',1.3);
hold on;
plot(Ang,v2t(:,3)','b','linewidth',1.3);
hold on;
plot(ang,v2(:,1)','r--','linewidth',1.3);
hold on;
plot(ang,v2(:,2)','g--','linewidth',1.3);
hold on;
plot(ang,v2(:,3)','b--','linewidth',1.3);

xlabel('Incidence angle','FontSize',12);
ylabel('{\itv_{p}}(m/s)','FontSize',12);
legend('{\itf} = 0.001 Hz','{\itf} = 10 Hz','{\itf} = 800 Hz');
xlim([0,360]);
set(gca,'XTick',[0 60 120 180 240 300 360]);


figure(4)
plot(Ang,inQ2(:,1),'r','linewidth',1.3);
hold on;
plot(Ang,inQ2(:,2),'g','linewidth',1.3);
hold on;
plot(Ang,inQ2(:,3),'b','linewidth',1.3);

xlabel('Incidence angle','FontSize',12);
ylabel('1/{\itQ_{p}}','FontSize',12);
legend('{\itf} = 0.1 Hz','{\itf} = 80 Hz','{\itf} = 2500 Hz');
xlim([0,360]);
set(gca,'XTick',[0 60 120 180 240 300 360]);

figure(5)
plot(Ang,v3t(:,1)','r','linewidth',1.3);
hold on;
plot(Ang,v3t(:,2)','g','linewidth',1.3);
hold on;
plot(Ang,v3t(:,3)','b','linewidth',1.3);
hold on;
plot(ang,v3(:,1)','r--','linewidth',1.3);
hold on;
plot(ang,v3(:,2)','g--','linewidth',1.3);
hold on;
plot(ang,v3(:,3)','b--','linewidth',1.3);

xlabel('Incidence angle','FontSize',12);
ylabel('{\itv_{p}}(m/s)','FontSize',12);
legend('{\itf} = 0.001 Hz','{\itf} = 10 Hz','{\itf} = 800 Hz');
xlim([0,360]);
set(gca,'XTick',[0 60 120 180 240 300 360]);


figure(6)
plot(Ang,inQ3(:,1),'r','linewidth',1.3);
hold on;
plot(Ang,inQ3(:,2),'g','linewidth',1.3);
hold on;
plot(Ang,inQ3(:,3),'b','linewidth',1.3);

xlabel('Incidence angle','FontSize',12);
ylabel('1/{\itQ_{p}}','FontSize',12);
legend('{\itf} = 0.1 Hz','{\itf} = 80 Hz','{\itf} = 2500 Hz');
xlim([0,360]);
set(gca,'XTick',[0 60 120 180 240 300 360]);