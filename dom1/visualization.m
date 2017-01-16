%visualization
figure;
hold on;
for i=0:5
in=['input_',num2str(i)];
out=['output_',num2str(i)];
loc=load(in);
val=load(out);
subplot(3,2,i+1);
cd ../util
plotSensors(loc,val);
cd ../dom1
end

%plot(test3_loc(:,1),test3_loc(:,2));

