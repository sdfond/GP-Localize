cd run4
tname = 'work'
fname = '1_2'
setname = 'D1_T'

O = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_OPITC', '_err'));
I = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_IOPITC', '_err'));
T = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_trunc', '_err'));
S = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_SoD', '_err'));
D = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_dr', '_err'));

step = length(O);
cd ..
so = sum(O)/step;
si = sum(I)/step;
ss = sum(S)/step;
st = sum(T)/step;
sd = sum(D)/step;


s = [st ss so si]
s = sd - s

bar(s)
set(gca,'XTickLabel',{'truncation', 'SoD', 'online PITC', 'improved PITC'})
ylabel('error');