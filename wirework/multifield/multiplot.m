%f=figure;
casenum = 3;
%hold on
A = zeros(421, 1);
B = zeros(421, 1);
C = zeros(421, 1);
D = zeros(421, 1);
E = zeros(421, 1);
F = zeros(421, 1);
G = zeros(421, 1);
for i=1:3
    str = strcat('run', num2str(i));
    cd(str)
    algo = 'OPITC';
    traj = 'D1_Twork';
    
    str = strcat(traj, '_F6_0_1_2_3_4_5_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    A = A + Z;
    
    str = strcat(traj, '_F1_0_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    B = B + Z;
    
    str = strcat(traj, '_F1_1_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    C = C + Z;
    
    str = strcat(traj, '_F1_2_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    D = D + Z;
    
    str = strcat(traj, '_F1_3_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    E = E + Z;
    
    str = strcat(traj, '_F1_4_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    F = F + Z;
    
    str = strcat(traj, '_F1_5_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    G = G + Z;
    cd ..
end
A = A / casenum;
B = B / casenum;
C = C / casenum;
D = D / casenum;
E = E / casenum;
F = F / casenum;
G = G / casenum;

% plot(A, 'k', 'LineWidth', 2);
% plot(B, 'g');
% plot(C, 'y');
% plot(D, 'b');
% plot(E, 'r');
% plot(F, 'c');
% plot(G, 'm');
% legend('Combine', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6');

st=1;
en = length(A);
%gd=ceil((en-st)/5);
gd=1;
x = [st:gd:en];
y = [A(x)'; B(x)'; C(x)'; D(x)'; E(x)'; F(x)'; G(x)'];
xl='Filtering/time step';
yl='Localization error';
%marker={'r*-','gs-','yv-','b^-', 'ks-', 'cx-', 'mo-'};
marker={'r-','g-','y-','b-', 'k-', 'c-', 'm-'};
legend={'Multiple fields', 'Field 1', 'Field 2', 'Field 3', 'Field 4', 'Field 5', 'Field 6'};
c=mFig(x,y,xl,yl,marker,legend);
c.lpos='NorthWest';
c.xlm=[0,430];
c.ylm=[0,45];

mPlot('plot',c,'wifi-multi');