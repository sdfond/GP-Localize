cd run1
O = load('real_OPITC_err');
I = load('real_IOPITC_err');
T = load('real_trunc_err');
S = load('real_SoD_err');
D = load('real_dr_err');
cd ..

step = length(O);
sum(D)/step
sum(T)/step
sum(S)/step
sum(O)/step
sum(I)/step



hold on
x = [1:step];
plot(x, T, '-k', 'MarkerSize', 10);
plot(x, S, '-b', 'MarkerSize', 10);
plot(x, O, '-r', 'MarkerSize', 10);
plot(x, I, '-m', 'MarkerSize', 10);

xlabel('step number');
ylabel('error');
legend('truncation', 'SoD', 'PITC', 'improved PITC', 'Location', 'NorthWest');
