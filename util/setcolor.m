%colorbar
function [cbar]=setcolor(cnum,ymin,ymax)
min_int=ceil(ymin*10)/10;
max_int=floor(ymax*10)/10;
%cc=jet(gd);

gd=ceil((max_int-min_int)/(cnum-2));
%gdt=ceil(180*100/(cnum-2))/100;
mi=27;
ma=193;
gdt=27.6;
% colorbar('YTickLabel',...
%     {'Freezing','Cold','Cool','Neutral',...
%      'Warm','Hot','Burning','Nuclear'});
%tick=[];
colormap(jet(cnum));
for i=1:cnum+1
    if i==1 
        cbar{i}='';
        tick(i)=0;
    elseif i==cnum+1
        cbar{i}='';
        tick(i)= mi+(i-2)*gdt;
    else
        cbar{i}=num2str(min_int+(i-2)*gd);
        tick(i)= mi+(i-2)*gdt;
    end
end
tick

% hcb = colorbar('YTickLabel',...
% {'Freezing','Cold','Cool','Neutral',...
% 'Warm','Hot','Burning','Nuclear'});
% set(hcb,'YTickMode','manual')
hcb=colorbar;
%set(hcb,'YTickMode','manual')
set(hcb,'YTickLabel',cbar);
%set(hcb,'YTick',tick);

% figure; 
% hold on;
% for i=1:8
%     plot([0 1],[0 i],'color',cc(i,:));
% end