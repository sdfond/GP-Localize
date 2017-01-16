run= 'run1';
file=[run,'/D1_Twork_F1_3_A4_0_1_2_3_OPITC_err'];
A = load(file);
x=1:100;
y=A(x);
xl='Step number';
yl='Error';
figure;
plot(x,y);
hold on;
x = [20:32];
plot(x, A(x), 'k-', 'LineWidth', 4);
x = [41:47];
plot(x, A(x), 'k--', 'LineWidth', 4);
x = [60:70];
plot(x, A(x), 'k:', 'LineWidth', 4);
ylim([0,10]);
xlabel(xl,'FontName','Times','FontSize', 20);
ylabel(yl,'FontName','Times','FontSize', 20);
set(gca,'FontSize',20);

cd ..
figure;
dispTrajOnMap(1, 3, 'work', 1, 100,'b','-',1);
dispTrajOnMap(1, 3, 'work', 20, 32,'k','-',4);
dispTrajOnMap(1, 3, 'work', 40, 47,'k','--',4);
dispTrajOnMap(1, 3, 'work', 60, 70,'k',':',4);
xlim([0,25]);
set(gca,'FontSize',20);

%plot(A(1:100));
% hold on
% x = [20:32];
% plot(x, A(x), 'k', 'LineWidth', 3);
% x = [41:47];
% plot(x, A(x), 'k', 'LineWidth', 3);
% x = [60:70];
% plot(x, A(x), 'k', 'LineWidth', 3);
% cd ..
