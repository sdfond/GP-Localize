casenum = 5
res = zeros(3, 3)
for i = [1,3,5,8,9]
    fdir = strcat('env1-', num2str(i));
    cd(fdir)
    for k = 0:1
        tname = 'working'
        fname = strcat('1_', num2str(k));
        setname = 'D0_T'
        O = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_OPITC', '_err'));
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
x = 1:3
y = res;
xl = ' ';
yl='Error';
marker={'-ko','-bs','-mv'};
legend={'Truncation', 'SoD',  'Improved online PITC'};

c=mFig(x,y,xl,yl,marker,legend);
c.baseline=0;
c.lpos='NorthWest';
c.xtk={'Temperature','Lighting','Humidity'};
c.ylm=[0,20]

mPlot('barx',c,'avgerror');