cd test
tname = 'real'
fname = '1_1'
setname = 'D2_T'

O = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_OPITC', '_err'));
I = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_IOPITC', '_err'));
T = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_trunc', '_err'));
S = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_SoD', '_err'));
D = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_dr', '_err'));
cd ..
hold on
x = [1:length(O)];
plot(x, T, '-k', 'MarkerSize', 10);
plot(x, S, '-b', 'MarkerSize', 10);
plot(x, O, '-r', 'MarkerSize', 10);
plot(x, I, '-m', 'MarkerSize', 10);


xlabel('step number');
ylabel('error');
legend('truncation', 'SoD', 'online FITC', 'improved PITC', 'Location', 'NorthWest');
