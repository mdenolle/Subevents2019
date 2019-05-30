function plot_simulations(xdat,tdat,pars,I_figure,FILE_name)

dt=tdat.Time(2)-tdat.Time(1);


[T_extend,X_extend]=get_rupture_extends(xdat);

STF=mean(xdat.SlipRate,1)*pars.L*pars.MU;


taus = pars.SIG0*pars.FRIC.MUs;
taud = pars.SIG0*pars.FRIC.MUd;

%% plot the results
close all
disp('Plotting ...')
figure(1)



% space-time plot of slip rate
subplot(5,5,[1 13])
%contourf(xdat.Time,xdat.X/1e3,xdat.SlipRate,10,'LineStyle','none')
imagesc(xdat.Time,xdat.X/1e3,xdat.SlipRate);
set(gca,'YDir','normal');
clb=colorbar('E'); % place colorbar INSIDE the contour plot to avoid axis
clb.Color='w';
% mis-aligned with the other subplots
title('Slip Rate (m/s)')
ylabel('X (km)')
xlim(T_extend)
ylim(X_extend/1e3)
% show also, for reference, the S wave front emanating from
% the hypocenter and a front 3 times slower


% moment rate
subplot(5,5,[16 18 21 23])
plot(xdat.Time,STF)
ylabel('Moment rate per unit width')
ylim([0 1.2]*max(STF))
xlim(T_extend)
xlabel('Time (s)')

hold on



%         % maximum slip rate, timeseries
%         subplot(5,5,[21 23])
%         plot(tdat.Time,tdat.MaxSlipRate);
%         xlabel('Time (s)')
%         ylabel('Max slip rate')

% moment rate spectral

T_long=-10:dt:1000;
MR_long=interp1(tdat.Time,tdat.MeanSlipRate*pars.L,T_long);
MR_long(isnan(MR_long))=0;
MR_long=MR_long/sum(MR_long*dt);
[MR_S,MR_f,~,~]=FFT_seimograph(MR_long,dt);

[fc_best,decay_best]=one_fc_fitting(MR_f(MR_f<1)',abs(MR_S(MR_f<1)),1);
MR_model=one_fc_model(MR_f,1,fc_best,decay_best);

subplot(5,5,[19 25])
loglog(MR_f,abs(MR_S),'-k');
hold on
loglog(MR_f(MR_f<1),MR_model(MR_f<1),'-r');
loglog(MR_f,MR_model,'--r');
loglog(fc_best,one_fc_model(fc_best,1,fc_best,decay_best),...
    '-ro','MarkerFaceColor','y','MarkerSize',15);

legend('Model','Brune-fit')
xlabel('Freq (Hz)')
grid on
xlim([10^(-4) ceil(1/dt)])
ylim([10^(-10) 10])

text(1e-3,1e-7,['f_c = ' num2str(fc_best,3) 'Hz'],'FontSize',15)
text(1e-3,1e-8,['n = ' num2str(decay_best,3)],'FontSize',15)

% slip snapshots
subplot(5,5,[4 14])


II_T=find(xdat.Time>T_extend(1) & xdat.Time<T_extend(2));
cmp_t=jet(length(II_T)+5);

for iit=1:round(length(II_T)/10):length(II_T)
    plot(xdat.Slip(:,II_T(iit)), xdat.X/1e3,'Color',cmp_t(II_T(iit),:));
    hold on
end
title('Slip (m)')
ylim(X_extend/1e3)

% max slip rate, spatial distribution
subplot(5,5,[5 15])
%plot(xdat.Stress(:,1)/1e6, xdat.X/1e3,'LineWidth',0.5)
Afill=fill(xdat.Stress(:,1)/1e6, xdat.X/1e3,'r');
Afill.LineStyle='none';
Afill.FaceAlpha=0.4;

xkm = [-1 1]*pars.L/2/1e3; % [xmin xmax] in km
hold on
plot([1 1]*pars.SIG0*pars.FRIC.MUs/1e6,xkm)
plot([1 1]*pars.SIG0*pars.FRIC.MUd/1e6,xkm)
title(['avg \Delta\tau=' num2str(pars.avg_stressdrop,2) 'MPa'])
xlim([taud/1e6-10 taus/1e6+10])
ylim(X_extend/1e3)

f1=figure(I_figure);
f1.Position=[50 50 1363 985];
f1.PaperUnits='points';
f1.PaperSize=f1.Position(3:4);
%pause
print('-dpng',FILE_name);