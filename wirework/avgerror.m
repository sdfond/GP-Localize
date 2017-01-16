casenum = 4
res = zeros(3, 6)
for i = 1:casenum
    fdir = strcat('run', num2str(i));
    cd(fdir)
    for k = 0:5
        tname = 'work';
        fname = strcat('1_', num2str(k));
        setname = 'D1_T';
        %O = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_OPITC', '_err'));
        I = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_IOPITC', '_err'));
        T = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_trunc', '_err'));
        S = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_SoD', '_err'));
        D = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_dr', '_err'));
        %res(:,k+1) = res(:,k+1) + [sum(D-T)/length(T); sum(D-S)/length(S); sum(D-I)/length(I)];
        res(:,k+1) = res(:,k+1) + [sum(T)/length(T); sum(S)/length(S); sum(I)/length(I)];
    end
    cd ..
end
res = res / casenum
%hold on
x = 1:6
y = res;
xl = ' ';
yl='Localization Error';
marker={'-ko','-bs','-mv'};
legend={'Truncation', 'SoD', 'GP-Localize'};

c=mFig(x,y,xl,yl,marker,legend);
c.baseline=0;
c.lpos='NorthWest';
c.xtk={'F1','F2','F3','F4','F5','F6'};
c.xlm=[0,7];
c.ylm=[0,35];


mPlot('barx',c,'avgerror');
