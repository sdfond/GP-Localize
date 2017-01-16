casenum = 1
res = zeros(4, 3)
for i = 5:5
    fdir = strcat('run', num2str(i));
    cd(fdir)
    for k = 0:2   
        tname = 'real'
        fname = strcat('1_', num2str(k));
        setname = 'D2_T'
        O = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_OPITC', '_err'));
        I = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_IOPITC', '_err'));
        T = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_trunc', '_err'));
        S = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_SoD', '_err'));
        D = load(strcat(setname, tname, '_F', fname, '_A4_0_1_2_3', '_dr', '_err'));
        %res(:,k+1) = res(:,k+1) + [sum(D-T)/length(T); sum(D-S)/length(S); sum(D-O)/length(O); sum(D-I)/length(I)];
        res(:,k+1) = res(:,k+1) + [sum(T)/length(T); sum(S)/length(S); sum(O)/length(O); sum(I)/length(I)];
    end
    cd ..
end
res = res / casenum
% hold on
% set(gca,'XTickLabel',str2mat('temperature', 'lighting', 'humidity'))
% plot(res(1,:), '-co');
% plot(res(2,:), '-bo');
% plot(res(3,:), '-r^');
% plot(res(4,:), '-m^');
% legend('truncation', 'SoD', 'online PITC', 'improved PITC');
% ylabel('improvedment');

x = 1:3
y = res;
xl = ' ';
yl='Improvement';
marker={'-ko','-bs','-r^','-mv'};
legend={'Truncation     ', 'SoD            ', 'Online PITC    ', 'Improved PITC  '};

c=mFig(x,y,xl,yl,marker,legend);

c.lpos='South';
c.baseline=-12;
c.xtk={'Temperature','Light','Humidity'};

mPlot('barx',c,'avgerror');

