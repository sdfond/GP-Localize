imagesc(layout);
hold on;

r=3.14/72;
rm=[cos(r) -1*sin(r);sin(r) cos(r)];
or=in*0+100;

plotSensors((in-or)*rm+or,out2);

