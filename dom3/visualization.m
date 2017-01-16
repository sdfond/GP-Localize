%visualization
figure;
hold on;
j=1;
for i=0:3:8
in=['input_',num2str(i)];
out=['output_',num2str(i)];
loc=load(in);
val=load(out);
subplot(2,2,j);
j=j+1;
cd ../util
plotSensors(loc,10*val);
cd ../dom3
end


