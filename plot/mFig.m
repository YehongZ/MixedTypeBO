%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAST MODIFICATION: 19 MAR. 2013
% ARIK JIE CHEN @ MAPLECG GROUP.
% SCHOOL OF COMPUTING, NUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize a figure
%x - vector/vectors of x
%y - vector/vectors of y
%leg - legend
%mker - marker
%xlb - label of x-axis
%ylb - label of y-axis
%lpos - position of legend
%ttl - title
%xlm - limit range to show for x-axis
%ylm - limit range to show for y-axis
%xtk - ticks to show in label of x-axis
%ytk - ticks to show in label of y-axis
function [c]=mFig(x,y,xl,yl,mark,leg)

sz=size(y);

if((sz(1)~=length(mark))|sz(1)~=length(leg))
disp('!!!ERROR: no. of markers and legends are different from no. of curves');
return;
end

if (sz(2)~=length(x))
disp('!!!ERROR: x and y have different points (columns)');
return;
end

%%%%%%%%%%%%%%%%%%%%%
c.x=x;
c.y=y;
c.leg=leg;
c.mker=mark;
c.xlb=xl;
c.ylb=yl;
%%%%%%%%%%%%%%%%%%%%%
c.lpos='';
c.ttl='';
c.xlm='';
c.ylm='';
c.xtk='';
c.ytk='';
