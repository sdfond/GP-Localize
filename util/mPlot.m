%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAST MODIFICATION: 19 MAR. 2013
% ARIK JIE CHEN @ MAPLECG GROUP.
% SCHOOL OF COMPUTING, NUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [p]=mPlot(type,dat,ofi)
% @desc Plot a figure
% @input:
%   type - plot/logx/logy/logxy/bar
%   dat  - a figure structure output by mFig
%   ofi  - file name for output file (empty string means no output file)
function [p]=mPlot(type,dat,ofi)


%%%%%predifine%%%%%%
msize=20; %marker size
lsize=20; %font size of legend
fsize=20; %font size of title/label
asize=20; %font size of axis ticks
fonts='Times'; %font
%%%%%%init%%%%%%%%%%%%%%
tt=dat.ttl;
lp=dat.lpos;
xlm=dat.xlm;
ylm=dat.ylm;
%%%%%%%%%%%%%%%%%%%%%%%
l=dat.leg;
lspec=dat.mker;
x=dat.x;
y=dat.y;
xl=dat.xlb;
yl=dat.ylb;
xtk=dat.xtk;
ytk=dat.ytk;
%baseline=dat.baseline;
%#curves=#rows
sz=size(y);
%if isempty(ofi)~=1
p=figure;
%end
box on;
hold on;

switch type    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %plot
  case 'plot'

    for j=1:sz(1)
        szx=size(x);
        szy=size(y);
        if szx(1)==1
            plot(x,y(j,:),lspec{j},'MarkerSize', msize);
        elseif szx(1)==szy(1)
            plot(x(j,:),y(j,:),lspec{j},'MarkerSize', msize);
        else
            disp('x is not compitable with y\n');
        end
    end
    set(gca,'FontSize',asize);
    %xlabel(xl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    %ylabel(yl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    xlabel(xl,'FontName',fonts,'FontSize', fsize);
    ylabel(yl,'FontName',fonts,'FontSize', fsize);
    legend(l);
    hleg=legend('show');
    if isempty(lp)==1
        lp='North';
    end
    
    if isempty(tt)~=1
        title(tt);
    end

    if isempty(xlm)~=1
        xlim(xlm);
    end

    if isempty(ylm)~=1
        ylim(ylm);
    end

    if isempty(xtk)~=1
        set(gca,'XTick',xtk);
    end

    if isempty(ytk)~=1
        set(gca,'YTick',ytk);
    end
    
    if isempty(lp)~=1
        set(hleg,'FontName',fonts,'FontSize',lsize,'Location',lp);
        set(hleg,'Interpreter','latex');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %semiloglog
   case 'logxy'
    hold off;
    for j=1:sz(1)
        szx=size(x);
        szy=size(y);
        if szx(1)==1
            loglog(x,y(j,:),lspec{j},'MarkerSize', msize);
        elseif szx(1)==szy(1)
            loglog(x(j,:),y(j,:),lspec{j},'MarkerSize', msize);
        else
            disp('x is not compitable with y\n');
        end
        if j==1
            hold on;
        end
    end
    set(gca,'FontSize',asize);
    %xlabel(xl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    %ylabel(yl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    xlabel(xl,'FontName',fonts,'FontSize', fsize);
    ylabel(yl,'FontName',fonts,'FontSize', fsize);
    legend(l);
    hleg=legend('show');
    if isempty(lp)==1
        lp='North';
    end
    
    if isempty(tt)~=1
        title(tt);
    end

    if isempty(xlm)~=1
        xlim(xlm);
    end

    if isempty(ylm)~=1
        ylim(ylm);
    end

    if isempty(xtk)~=1
        set(gca,'XTick',xtk);
    end

    if isempty(ytk)~=1
        set(gca,'YTick',ytk);
    end
    
    if isempty(lp)~=1
        set(hleg,'FontName',fonts,'FontSize',lsize,'Location',lp);
        set(hleg,'Interpreter','latex');
    end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %semilogx
   case 'logx'
    hold off;
    for j=1:sz(1)
        szx=size(x);
        szy=size(y);
        if szx(1)==1
            semilogx(x,y(j,:),lspec{j},'MarkerSize', msize);
        elseif szx(1)==szy(1)
            semilogx(x(j,:),y(j,:),lspec{j},'MarkerSize', msize);
        else
            disp('x is not compitable with y\n');
        end
        if j==1
            hold on;
        end
    end
    set(gca,'FontSize',asize);
    %xlabel(xl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    %ylabel(yl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    xlabel(xl,'FontName',fonts,'FontSize', fsize);
    ylabel(yl,'FontName',fonts,'FontSize', fsize);
    legend(l);
    hleg=legend('show');
    if isempty(lp)==1
        lp='North';
    end
    
    if isempty(tt)~=1
        title(tt);
    end

    if isempty(xlm)~=1
        xlim(xlm);
    end

    if isempty(ylm)~=1
        ylim(ylm);
    end

    if isempty(lp)~=1
        set(hleg,'FontName',fonts,'FontSize',lsize,'Location',lp);
        set(hleg,'Interpreter','latex');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%semilogy
   case 'logy'
    hold off;
    for j=1:sz(1)
        szx=size(x);
        szy=size(y);
        if szx(1)==1
            semilogy(x,y(j,:),lspec{j},'MarkerSize', msize);
        elseif szx(1)==szy(1)
            semilogy(x(j,:),y(j,:),lspec{j},'MarkerSize', msize);
        else
            disp('x is not compitable with y\n');
        end
        if j==1
            hold on;
        end
    end
    set(gca,'FontSize',asize);
    %xlabel(xl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    %ylabel(yl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    xlabel(xl,'FontName',fonts,'FontSize', fsize);
    ylabel(yl,'FontName',fonts,'FontSize', fsize);
    legend(l);
    hleg=legend('show');
    if isempty(lp)==1
        lp='North';
    end
    
    if isempty(tt)~=1
        title(tt);
    end

    if isempty(xlm)~=1
        xlim(xlm);
    end

    if isempty(ylm)~=1
        ylim(ylm);
    end

    if isempty(xtk)~=1
        set(gca,'XTick',xtk);
    end

    if isempty(ytk)~=1
        set(gca,'YTick',ytk);
    end

    
    if isempty(lp)~=1
        set(hleg,'FontName',fonts,'FontSize',lsize,'Location',lp);
        set(hleg,'Interpreter','latex');
    end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %bar
  case 'bar'
    
    bar(y');
    set(gca,'FontSize',asize);
    set(gca, 'XTick', floor(x));
    %settick('x',floor(x),fonts,fsize);
    %xlabel(xl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    %ylabel(yl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    xlabel(xl,'FontName',fonts,'FontSize', fsize);
    ylabel(yl,'FontName',fonts,'FontSize', fsize);
    colormap summer;
    
    if isempty(tt)~=1
        title(tt);
    end

    %xlim(xlm);
    if isempty(l)~=1
      legend(l);
      hleg=legend('show');
      if isempty(lp)
        lp='North';
      end
      set(hleg,'FontName',fonts,'FontSize',lsize,'Location',lp);
      set(hleg,'Interpreter','latex');
    end

    case 'barx'
    
    bar(y');
    set(gca,'FontSize',asize);
    set(gca, 'XTick', floor(x));
    settick('x',xtk, dat.baseline,fsize);
    %xlabel(xl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    %ylabel(yl,'FontName',fonts,'FontSize', fsize,'Interpreter','latex');
    xlabel(xl,'FontName',fonts,'FontSize', fsize);
    ylabel(yl,'FontName',fonts,'FontSize', fsize);
    colormap summer;
    
    if isempty(tt)~=1
        title(tt);
    end

    
    if isempty(xlm)~=1
        xlim(xlm);
    end

    if isempty(ylm)~=1
        ylim(ylm);
    end
    
    %xlim(xlm);
    if isempty(l)~=1
      legend(l);
      hleg=legend('show');
      if isempty(lp)
        lp='North';
      end
      set(hleg,'FontName',fonts,'FontSize',lsize,'Location',lp);
      set(hleg,'Interpreter','latex');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    disp(['bad type: try plot|bar']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(ofi)~=1
%ofi=[of,num2str(i)];
print('-depsc',ofi);
close(p);

end
end

%set(haxes(1),'XTick',[1 10 100 ...])
