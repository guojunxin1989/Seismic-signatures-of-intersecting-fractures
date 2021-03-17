clear;
close;
clc;

load parallel_case f v inQ;
f1=f;
v1=v;
inQ1=inQ;

load parallel_case_Guo_2020 f v inQ;
f2=f;
v2=v;
inQ2=inQ;

figure(1)
semilogx(f1,v1(:,1),'m','linewidth',1.3);
hold on;
semilogx(f1,v1(:,2),'r','linewidth',1.3);
hold on;
semilogx(f1,v1(:,3),'g','linewidth',1.3);
hold on;
semilogx(f1,v1(:,4),'b','linewidth',1.3);
hold on;
semilogx(f1,v1(:,5),'k','linewidth',1.3);
hold on;
semilogx(f2,v2(:,1),'mo','MarkerSize',5);
hold on;
semilogx(f2,v2(:,2),'ro','MarkerSize',5);
hold on;
semilogx(f2,v2(:,3),'go','MarkerSize',5);
hold on;
semilogx(f2,v2(:,4),'bo','MarkerSize',5);
hold on;
semilogx(f2,v2(:,5),'ko','MarkerSize',5);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('{\itP}-wave velocity (m/s)','FontSize',12);
legend('\theta = 0','\theta = 30','\theta = 45','\theta = 60','\theta = 90');
set(gca, 'FontSize', 12);
xlim([0.001 100000]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);

figure(2)
loglog(f1,inQ1(:,1),'m','linewidth',1.3);
hold on;
loglog(f1,inQ1(:,2),'r','linewidth',1.3);
hold on;
loglog(f1,inQ1(:,3),'g','linewidth',1.3);
hold on;
loglog(f1,inQ1(:,4),'b','linewidth',1.3);
hold on;
loglog(f1,inQ1(:,5),'k','linewidth',1.3);
hold on;
loglog(f2,inQ2(:,1),'mo','MarkerSize',5);
hold on;
loglog(f2,inQ2(:,2),'ro','MarkerSize',5);
hold on;
loglog(f2,inQ2(:,3),'go','MarkerSize',5);
hold on;
loglog(f2,inQ2(:,4),'bo','MarkerSize',5);
hold on;
loglog(f2,inQ2(:,5),'ko','MarkerSize',5);

xlabel('Frequency (Hz)','FontSize',12);
ylabel('Q_{p}^{-1}','FontSize',12);
legend('\theta = 0','\theta = 30','\theta = 45','\theta = 60','\theta = 90');
set(gca, 'FontSize', 12);
xlim([0.001 100000]);
ylim([0.000001 0.1]);
set(gca,'XTick',[0.001 0.01 0.1 1 10 100 1000 10000 100000]);

