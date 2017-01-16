loc = load('real_loc');
val = load('real_obs');
cd ../util/
plotSensors(loc(:,1:2), val(:,2)*10)
cd ../traj2
hold on
plot(loc(180:210,1), loc(180:210,2), '-k')