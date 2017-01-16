%f=figure;
casenum = 5;
%hold on
%A = zeros(421, 1);
B = zeros(336, 1);
C = zeros(336, 1);
%D = zeros(421, 1);
E = zeros(336, 1);
%F = zeros(421, 1);
%G = zeros(421, 1);
for i=[1,3,5,8,9]
    str = strcat('env1-', num2str(i));
    cd(str)
    algo = 'IOPITC';
    traj = 'D0_Tworking';
    
    str = strcat(traj, '_F1_0_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    B = B + Z;
    
    str = strcat(traj, '_F1_1_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    C = C + Z;
    
    cd ..
end

for  i =[4,7:9]
    str = strcat('env1-', num2str(i));
    cd(str)
    algo = 'IOPITC';
    traj = 'D0_Tworking';
    
    str = strcat(traj, '_F2_0_1_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    E = E + Z;
    cd ..
end
%A = A / casenum;
B = B / casenum;
C = C / casenum;
%D = D / casenum;
E = E / 4;
%F = F / casenum;
%G = G / casenum;
sum(B)/length(B)
sum(C)/length(C)
sum(E)/length(E)
% plot(A, 'k', 'LineWidth', 2);
% plot(B, 'g');
% plot(C, 'y');
% plot(D, 'b');
% plot(E, 'r');
% plot(F, 'c');
% plot(G, 'm');
% legend('Combine', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6');

st=1;
en = length(E);
%gd=ceil((en-st)/16);
gd=1;
x = [st:gd:en];
%y = [A(x)'; B(x)'; C(x)'; D(x)'; E(x)'; F(x)'; G(x)'];
y = [E(x)'; B(x)'; C(x)'];
xl='Filtering/time step';
yl='Localization error';
%marker={'r*-','gs-','b^-'};
marker={'r', 'b', 'g'};
%marker={'r-','g-','y-','b-', 'k-', 'c-', 'm-'};
%legend={'Combine', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6'};
legend={'Multiple fields', 'Temperature field', 'Light field'};
c=mFig(x,y,xl,yl,marker,legend);
c.lpos='NorthWest';
c.xlm=[0,350];
c.ylm=[0,7.8];

mPlot('plot',c,'env-multi');