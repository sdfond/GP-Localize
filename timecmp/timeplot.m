T = load('D2_Treal_F1_0_A6_0_1_2_3_4_5_trunc_time');
S = load('D2_Treal_F1_0_A6_0_1_2_3_4_5_SoD_time');
F = load('D2_Treal_F1_0_A6_0_1_2_3_4_5_FGP_time');
I = load('D2_Treal_F1_0_A6_0_1_2_3_4_5_OPITC_time');
P = load('D2_Treal_F1_0_A6_0_1_2_3_4_5_PITC_time');
O = load('D2_Treal_F1_0_A6_0_1_2_3_4_5_IOPITC_time');
T = T / 1000;
S = S / 1000;
F = F / 1000;
O = O / 1000;
P = P / 1000;
I = I / 1000;
% T = log(T);
% S = log(S);
% F = log(F);
% I = log(I);
% O = log(O);
% P = log(P);
%hold on
% plot(T, 'g');
% plot(S, 'c');
% plot(O, 'r');
% plot(I, 'm');
% plot(F, 'k');
% plot(P, 'b');
st=1;
en=length(T);
gd=ceil((en-st)/10);
x=[st:gd:en];
y=[T(x)'; S(x)';  F(x)'; P(x)'; I(x)'];
xl='Filtering/time steps';
yl='Incurred time (s) per filtering step';
marker={'gx-','ks-','bv-','m^-', 'ro-'};
legend={'SoD-Truncate', 'SoD-Even', 'Full GP', 'Offline PITC',  'GP-Localize'};
c=mFig(x,y,xl,yl,marker,legend);
c.lpos='NorthWest';

mPlot('logy',c,'timeplot');

% xlabel('step number');
% ylabel('running time (log second)');
%legend('trunation', 'SoD', 'online PITC', 'improved online PITC', 'FGP', 'batch PITC');