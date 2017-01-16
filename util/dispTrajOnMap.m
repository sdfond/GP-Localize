%visualization
function [p]=dispTrajOnMap(dataset,field,tname,st,en,lc,ls,lw)
%figure;
hold on;
%dataset=1;
dom=['dom',num2str(dataset)];
trj=['traj',num2str(dataset)];
%field=1;
%tname='test1';
%for i=0:5
in=[dom,'/input_',num2str(field)];
out=[dom,'/output_',num2str(field)];
loc=load(in);
val=load(out);
%subplot(3,2,i+1);
%cd ../util
plotSensors(loc,val);
%cd ../dom1
%end
tf=[trj,'/',tname,'_loc'];
tloc=load(tf);
l=[lc,ls];
p=plot(tloc(st:en,1),tloc(st:en,2),l,'LineWidth',lw);
l=[lc,'x'];
plot(tloc(st,1),tloc(st,2),l,'MarkerSize',20,'LineWidth',lw);
l=[lc,'o'];
plot(tloc(en,1),tloc(en,2),l,'MarkerSize',20,'LineWidth',lw);