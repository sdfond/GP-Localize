cd run2
tname = 'work'
fname = '1_1'
setname = 'D1_T'

O = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_OPITC', '_err'));
I = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_IOPITC', '_err'));
T = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_trunc', '_err'));
S = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_SoD', '_err'));
D = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_dr', '_err'));
cd ..
%hold on
st=1;
en = length(O);
gd=ceil((en-st)/10);
x = [st:gd:en];
y = [T(x)'; S(x)'; O(x)'; I(x)'];
xl='Step number';
yl='Error';
marker={'-ko','-bs','-r^','-mv'};
legend={'Truncation', 'SoD', 'Online PITC', 'Improved PITC'};
c=mFig(x,y,xl,yl,marker,legend);
c.lpos='NorthWest';

mPlot('plot',c,'');

% plot(x, T, '-k', 'MarkerSize', 10);
% plot(x, S, '-b', 'MarkerSize', 10);
% plot(x, O, '-r', 'MarkerSize', 10);
% plot(x, I, '-m', 'MarkerSize', 10);
% 
% 
% xlabel('step number');
% ylabel('error');
% legend('truncation', 'SoD', 'online FITC', 'improved PITC', 'Location', 'NorthWest');
