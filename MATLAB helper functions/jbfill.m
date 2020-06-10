function[fillhandle,msg]=jbfill(xpoints,upper,lower,color,edge,add,transparency,stair)
%USAGE: [fillhandle,msg]=jbfill(xpoints,upper,lower,color,edge,add,transparency,stair)
%This function will fill a region with a color between the two vectors provided
%using the Matlab fill command.
%
%fillhandle is the returned handle to the filled region in the plot.
%xpoints= The horizontal data points (ie frequencies). Note length(Upper)
%         must equal Length(lower)and must equal length(xpoints)!
%upper = the upper curve values (data can be less than lower)
%lower = the lower curve values (data can be more than upper)
%color = the color of the filled area 
%edge  = the color around the edge of the filled area
%add   = a flag to add to the current plot or make a new one.
%transparency is a value ranging from 1 for opaque to 0 for invisible for
%the filled color only.
%
%John A. Bockstege November 2006;
%Example:
%     a=rand(1,20);%Vector of random data
%     b=a+2*rand(1,20);%2nd vector of data points;
%     x=1:20;%horizontal vector
%     [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,rand(1,1))
%     grid on
%     legend('Datr')
if nargin<8;stair=0;end %default is to have smooth curve, not stairs
if nargin<7;transparency=.5;end %default is to have a transparency of .5
if nargin<6;add=1;end     %default is to add to current plot
if nargin<5;edge='k';end  %dfault edge color is black
if nargin<4;color='b';end %default color is blue

%-----------------
% DA added
upper(find(isnan(upper)))=0;
lower(find(isnan(lower)))=0;
if edge == -1;
    edge = 'k';
    noedge = 1;
else
    noedge = 0;
end
if stair
    if diff(size(xpoints))<1 xpoints = xpoints'; end
    if diff(size(upper))<1 upper = upper'; end
    if diff(size(lower))<1 lower = lower'; end
    dx = mean(diff(xpoints));
    xpoints = reshape(addortho([-dx/2; dx/2],xpoints),[],1)';
    upper = reshape(repmat(upper,2,1),[],1)';
    lower = reshape(repmat(lower,2,1),[],1)';
end
%------------------

if length(upper)==length(lower) && length(lower)==length(xpoints)
    msg='';
    filled=[upper,fliplr(lower)];
    xpoints=[xpoints,fliplr(xpoints)];
    if add
        hold on
    end
    fillhandle=fill(xpoints,filled,color);%plot the data
    set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
    
    %-----------
    % DA added
    if noedge
        set(fillhandle,'LineStyle','none');
    end
    
    if add
        hold off
    end
    set(gca,'Layer','top');
else
    msg='Error: Must use the same number of points in each vector';
end
