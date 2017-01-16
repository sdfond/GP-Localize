%visualization
figure;
hold on;
for i=0:2
in=['input_',num2str(i)];
out=['output_',num2str(i)];
loc=load(in);
val=load(out);
subplot(3,1,i+1);
cd ../util
plotSensors(loc,val*10);
cd ../dom0
end

%plot(test3_loc(:,1),test3_loc(:,2));

