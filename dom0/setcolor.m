%colorbar
function [cbar]=setcolor(cnum,ymin,ymax)
min_int=ceil(ymin);
max_int=floor(ymax);
%cc=jet(gd);

gd=ceil((max_int-min_int)/(cnum-2))


% colorbar('YTickLabel',...
%     {'Freezing','Cold','Cool','Neutral',...
%      'Warm','Hot','Burning','Nuclear'});
 
colormap(jet(cnum));
for i=1:cnum+1
    if i==1 | i==cnum+1
        cbar{i}='';
    else
        cbar{i}=num2str(min_int+(i-2)*gd);
    end
end

% hcb = colorbar('YTickLabel',...
% {'Freezing','Cold','Cool','Neutral',...
% 'Warm','Hot','Burning','Nuclear'});
% set(hcb,'YTickMode','manual')

hcb = colorbar('YTickLabel',cbar);

% figure; 
% hold on;
% for i=1:8
%     plot([0 1],[0 i],'color',cc(i,:));
% end