%visualization
load work_loc
st=1;
en=length(work_loc);
for i=0:5
hfig=figure;
hold on;
%cd ../dom0
in=['../dom1/input_',num2str(i)];
out=['../dom1/output_',num2str(i)];
loc=load(in);
val=load(out);


%subplot(3,1,i+1);
%cd ../util
plotSensors(loc,val);
% xlim([-30,30]);
% ylim([-6,6]);
plot(work_loc(st:en,1),work_loc(st:en,2),'k');
plot(work_loc(st,1),work_loc(st,2),'kx','MarkerSize',20,'LineWidth',1);
plot(work_loc(en,1),work_loc(en,2),'ko','MarkerSize',20,'LineWidth',1);

set(hfig, 'Position', [0 0 1060 240])
ylim([-6,6])
xlim([-24,24])
set(gca,'YTick',[-6:2:6])
set(gca,'XTick',[-24:2:24])
fname=['wifi-trj-f',num2str(i)];
set(hfig, 'paperpositionmode', 'auto');
print('-depsc',fname);
close(hfig);
end


%cd ../traj0
