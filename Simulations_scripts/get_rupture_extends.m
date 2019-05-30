function [T_extend,X_extend]=get_rupture_extends(xdat)
%% get the rupture extension from the model
% [T_extend,X_extend]=get_rupture_extends(tdat,xdat)

II_time=1:min((find(sum(xdat.SlipRate,1)>0, 1, 'last' )),length(sum(xdat.SlipRate,1)));
II_X=(find(xdat.Slip(:,end)>0, 1 )-5):(find(xdat.Slip(:,end)>0, 1, 'last' )+5);

Temporal=xdat.Time(II_time);
Spatial=xdat.X(II_X);

T_extend=[min(Temporal) max(Temporal)];
X_extend=[min(Spatial) max(Spatial)];