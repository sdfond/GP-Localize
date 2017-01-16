%visualization
load working_loc
st=1;
en=length(working_loc);
for i=1:1
figure;
hold on;
%cd ../dom0
in=['../dom0/input_',num2str(i)];
out=['../dom0/output_',num2str(i)];
loc=load(in);
val=load(out);

if i == 0
    val = val * 1.8 + 32;
end
%subplot(3,1,i+1);
%cd ../util
plotSensors(loc,val);

plot(working_loc(st:en,1),working_loc(st:en,2),'k');
plot(working_loc(st,1),working_loc(st,2),'kx','MarkerSize',30,'LineWidth',1);
plot(working_loc(en,1),working_loc(en,2),'ko','MarkerSize',30,'LineWidth',1);
end
xlim([-11,11])
ylim([-9,9])
set(gca,'YTick',[-8:2:8])
set(gca,'XTick',[-12:2:12])
set(hfig, 'Position', [0 0 560 360])
print('-depsc','env-light');
%cd ../traj0
