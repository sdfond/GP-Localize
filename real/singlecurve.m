cd run3
tname = 'real'
setname = 'D2_T'

hold on
O1 = load(strcat(setname, tname, '_F', '1_0', '_A4_0_1_2_3', '_OPITC', '_err'));
plot(O1, '-r');
O2 = load(strcat(setname, tname, '_F', '1_1', '_A4_0_1_2_3', '_OPITC', '_err'));
plot(O2, '-c');
O3 = load(strcat(setname, tname, '_F', '1_2', '_A4_0_1_2_3', '_OPITC', '_err'));
plot(O3, '-b');

legend('temperature', 'lighting', 'humidity');
cd ..
