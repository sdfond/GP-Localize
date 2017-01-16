function [pp]=plotSensors(in,out)
out=out;
colnum=6;
figure;
hold on;
min(out)
max(out)
for i=1:length(out)
c=getcolor(colnum,min(out),max(out),out(i));
plot(in(i,1),in(i,2),'o','MarkerFaceColor',c)
end
setcolor(colnum,min(out),max(out));
