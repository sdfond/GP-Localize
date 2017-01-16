function settick(axis,ticks,baseline, fsize)
%fsize=32;
n=length(ticks);
tkx=get(gca,'XTick');tky=get(gca,'YTick');
switch axis
    case 'x'
        w=linspace(tkx(1),tkx(end),n);
        set(gca, 'XTick', w, 'XTickLabel', []);
        %yh=(15*w(1)-w(end))/14
        %yh=baseline;
        for i=1:n
            text('Interpreter','tex','String',ticks{i},'Position',[w(i),baseline],'horizontalAlignment', 'center','verticalAlignment', 'top','FontName','Times','FontSize', fsize);
        end
    case 'y'
        w=linspace(tky(1),tky(end),n);
        set(gca, 'YTick', w, 'YTickLabel', []);
        xh=(11*w(1)-w(end))/10;
        for i=1:n
            text('Interpreter','tex','String',ticks{i},'Position',[xh,w(i)],'horizontalAlignment', 'center');
        end
end
