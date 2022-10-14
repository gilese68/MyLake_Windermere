function harmboxplot(med,q1,q3,upadj,loadj,lengthx,yy,notch,sym,vert,day)

%JFBOXPLOT(MED,Q1,Q3,UPADJ,LOADJ,LENGTHX,YY,NOTCH,SYM,VERT)
%
% MED is an array of the medians for each set of data.
% Q1 is the lower quartile for those sets.
% Q3 is the upper quartile for these sets.
%
% UPADJ is an array of same size, and for 
% each column of data is given by:
%vhi = q3+1.5*(q3-q1);
%upadj = max(x(x<=vhi));
%if (isempty(upadj)), upadj = q3; end
%
% Where whis is as for BOXPLOT (default is 1.5 
% if you wish to use it)
%
% LOADJ for each column given by:
% vlo = q1-whis*(q3-q1);
% loadj = min(x(x>=vlo));
% if (isempty(loadj)), loadj = q1; end
%
% LENGTHX is the length of each column - should be an
% array the same size as med
%
% YY is a cell array.  The entry for each column
% is given by:
% yy = x(x<loadj | x > upadj);
% If you don't care about this data, set each cell
% to be []
%
% NOTCH, SYM and VERT are as for BOXPLOT but must be specified.
% The defaults for BOXPLOT are 0, 'r+' and 1
% 
% MED, Q1, Q3, UPADJ, LOADJ, LENGTHX and YY must have the
% same dimension.  For the ith column of data, median
% is MED(i), upper quartile is Q1(i), etc.
%
%
% EXAMPLE:  
% %Prepare data:
% x=rand(100,6);
% whis=1.5;
% for i=1:6
%   xt=x(:,i);
%   med(i)= prctile(xt,50);
%   q1(i) = prctile(xt,25);
%   q3(i) = prctile(xt,75);
%   vhi(i) = q3(i)+whis*(q3(i)-q1(i));
%   upadj(i) = max(xt(xt<=vhi(i)));
%  if (isempty(upadj(i))), upadj(i) = q3(i); end
%   vlo(i) = q1(i)-whis*(q3(i)-q1(i));
%   loadj(i) = min(xt(xt>=vlo(i)));
%   if (isempty(loadj(i))), loadj(i) = q1(i); end
%   lengthx(i)=100;
%   yy{i} = xt(xt<loadj(i) | xt > upadj(i));
% end
%
% %Run file
% jfboxplot(med,q1,q3,upadj,loadj,lengthx,yy,1,'r+',1);

%cla
set(gca,'NextPlot','add','Box','on');

n=length(med);
if vert
    set(gca,'YLabel',text(0,0,'Values'));
    set(gca,'XLabel',text(0,0,' ')); 
else
    set(gca,'XLabel',text(0,0,'Values'));
    set(gca,'YLabel',text(0,0,'Site')); 
end

lf=n*min(0.15,0.5/n);

for i=1:n
    boxutilnew(med(i),q1(i),q3(i),upadj(i),loadj(i),lengthx(i),yy{i},notch,i,lf,sym,vert,day);
    
    % tp 21/oct/04
%      hold on
%      text(i-0.075,obse(SITES(i)),'o')
%      hold on
     %text(i-0.25,0,sitelab(SITES(i)),'rotation',90)
     
     %text((1:length(obse))-0.25,obse','o');
end


set(gca,'NextPlot','replace');


function boxutilnew(med,q1,q3,upadj,loadj,lengthx,yy,notch,lb,lf,sym,vert,day)

x1 = lb*ones(1,2);
x1(:) = day;
lb = day;
lf = 5;
x2 = x1+[-0.25*lf,0.25*lf];

if length(yy)==0
   yy = loadj;
   [a1 a2 a3 a4] = colstyle(sym);
   sym = [a2];
end

xx = lb*ones(1,length(yy));
    lbp = lb + 0.5*lf;
    lbm = lb - 0.5*lf;


upadj = max(upadj,q3);
loadj = min(loadj,q1);

% Set up (X,Y) data for notches if desired.
if ~notch
    xx2 = [lbm lbp lbp lbm lbm];
    yy2 = [q3 q3 q1 q1 q3];
    xx3 = [lbm lbp];
else
    n1 = med + 1.57*(q3-q1)/sqrt(lengthx);
    n2 = med - 1.57*(q3-q1)/sqrt(lengthx);
    if n1>q3, n1 = q3; end
    if n2<q1, n2 = q1; end
    lnm = lb-0.25*lf;
    lnp = lb+0.25*lf;
    xx2 = [lnm lbm lbm lbp lbp lnp lbp lbp lbm lbm lnm];
    yy2 = [med n1 q3 q3 n1 med n2 q1 q1 n2 med];
    xx3 = [lnm lnp];
end
yy3 = [med med];

% Determine if the boxes are vertical or horizontal.
% The difference is the choice of x and y in the plot command.
if vert
    plot(x1,[q3 upadj],'k-',x1,[loadj q1],'k-',...
        x2,[loadj loadj],'k-',...
        x2,[upadj upadj],'k-',xx2,yy2,'k-',xx3,yy3,'k-',xx,yy,sym)
else
    plot([q3 upadj],x1,'k-',[loadj q1],x1,'k-',...
        [loadj loadj],x2,'k-',...
        [upadj upadj],x2,'k-',yy2,xx2,'k-',yy3,xx3,'k-',yy,xx,sym)
end
