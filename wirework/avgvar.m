casenum = 4
res = zeros(4, 6)
for i = 1:casenum
    fdir = strcat('run', num2str(i));
    cd(fdir)
    for k = 0:5
        tname = 'work';
        fname = strcat('1_', num2str(k));
        setname = 'D1_T';
        %O = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_OPITC', '_var'));
        I = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_IOPITC', '_var'));
        T = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_trunc', '_var'));
        S = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_SoD', '_var'));
        D = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_dr', '_var'));
        %res(:,k+1) = res(:,k+1) + [sum(D-T)/length(T); sum(D-S)/length(S); sum(D-I)/length(I)];
        I = sqrt(I(:,1).^2+I(:,2).^2);
        T = sqrt(T(:,1).^2+T(:,2).^2);
        S = sqrt(S(:,1).^2+S(:,2).^2);
        D = sqrt(D(:,1).^2+D(:,2).^2);
        res(:,k+1) = res(:,k+1) + [sum(D)/length(D); sum(T)/length(T); sum(S)/length(S); sum(I)/length(I)];
    end
    cd ..
end
res = res / casenum
%hold on
%x = 1:6
%y = res;
%xl = ' ';
%yl='Localization varor';
%marker={'-ko','-bs','-mv'};
%legend={'Truncation', 'SoD', 'GP-Localize'};

%c=mFig(x,y,xl,yl,marker,legend);
%c.baseline=0;
%c.lpos='NorthWest';
%c.xtk={'F1','F2','F3','F4','F5','F6'};
%c.xlm=[0,7];
%c.ylm=[0,35];


mPlot('barx',c,'avgvar');
