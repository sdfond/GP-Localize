%1,4,5
f=figure;
casenum = 2;
hold on
A = zeros(561, 1);
B = zeros(561, 1);
C = zeros(561, 1);
D = zeros(561, 1);
E = zeros(561, 1);
F = zeros(561, 1);
G = zeros(561, 1);
K = zeros(561, 1);

for i=[4,5]
    str = strcat('run', num2str(i));
    cd(str)
    algo = 'IOPITC';
    traj = 'D2_Treal';
    
    str = strcat(traj, '_F3_0_1_2_A4_0_1_2_3_', algo, '_err');
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
    
    str = strcat(traj, '_F2_0_1_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    E = E + Z;
    
    str = strcat(traj, '_F2_1_2_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    F = F + Z;
    
    str = strcat(traj, '_F2_0_2_A4_0_1_2_3_', algo, '_err');
    Z = load(str);
    G = G + Z;
    
    str = strcat(traj, '_F2_0_1_A4_0_1_2_3_', 'dr', '_err');
    Z = load(str);
    K = K + Z;
    
    cd ..
end
A = A / casenum;
B = B / casenum;
C = C / casenum;
D = D / casenum;
E = E / casenum;
F = F / casenum;
G = G / casenum;
K = K / casenum;
sum(G)/length(G)
sum(K)/length(K)

plot(A, 'k', 'LineWidth', 2);
plot(B, 'r');
plot(C, 'b');
plot(D, 'c');
plot(E, 'g', 'LineWidth', 2);
plot(F, 'm', 'LineWidth', 2);
plot(G, 'y', 'LineWidth', 2);
legend('combine', '1', '2', '3', '1+2', '2+3', '1+3');