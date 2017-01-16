%%%%%%%%%%%%%%%% 1.Process the raw results into row vectors %%%%%%%%%%%%%%%
x=[1:5]; %  x axis
y=[1:5;2:6;3:7]; %  y axis
%%%%%%%%%%%%%%%% 2. initialize a figure structure %%%%%%%%%%%%%%%%%%%%%%%%%
xl='This is x label'; % 1. xlabel
yl='This is y label'; % 2. ylabel
% 3. legend the number of strings should be the 
%   same as the number of curves in yl
legend={'Legend $\alpha$','Legend $\beta$','Legend $\gamma$'}; 

% 4. marker is a string of three fields including
%   i)    line style (i.e., solid line '-', dash line '--')
%   ii)   color (i.e., 'm','r','g','b','k','y')
%   iii)  marker style (i.e., 'o','v','^','s','x',',')
% use 'doc plot' in Matlab for more information
marker={'-bo','-mv','-r^'};  

c=mFig(x,y,xl,yl,marker,legend);

%%%%%%%%%%%%%%% 2.5 Adjust legend/range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the legend: 'North', 'South', 'West', 'East', 'NorthWest' ...
c.lpos='NorthWest'; 
%c.ttl='No title'; % title
c.xlm=[0.5,6]; % range of plotted x 
c.ylm=[0,8]; % range of plotted y 
c.xtk=[0.5,1,2,4,6]; % ticks of x axis
c.ytk=[0,2,4,6,8]; % ticks of y axis

%%%%%%%%%%%%%% 3. plot/replot the figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the curve: plot/logx/logy/logxy/bar
mPlot('plot',c,'file');
% if the last argument is set empty, the figure is output to screen
mPlot('plot',c,'');
