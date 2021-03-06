%colorbar
function [c]=getcolor(cnum,ymin,ymax,y)

min_int=ceil(ymin);
max_int=floor(ymax);
cc=jet(cnum);

gd=(max_int-min_int)/(cnum-2)
ci=floor((y-min_int)/gd)+2

c=cc(ci,:);
% figure; 
% hold on;
% for i=1:8
%     plot([0 1],[0 i],'color',cc(i,:));
% end