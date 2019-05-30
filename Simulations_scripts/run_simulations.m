cd /Users/Yin9xun/Dropbox/Rupture_dynamic_codes/SBIEMLAB

clear
close all


Index_test=0;
plot_movie=0;
CorrL1=[10 5 2];
Dc_list=[0.4 0.8 1.2];



for i_Dc=1:length(Dc_list)
    %%
    
    rng(31) % for the random background stress
    Phs_BG=2*pi*rand(1,1024*2);
    
    % different random seeds to generate heterogeneous pre stress
    % rng(9) % for i_Phs=1:300 
    rng(13) % for i_Phs=301:600
    
    for i_Phs=301:600
        
        close all
        
        
        for iCorrL=10%CorrL1
            for iHurst=0.8 %0 corresponds to pre-setted stress (non-random)
                
                tic
                %pause(3)
                
                %close all
                clc
                %CorrL=40;              
                output_dir=['/Users/Yin9xun/Work/Two_Asperities_scaling/Dc_finer_' num2str(Dc_list(i_Dc))];
                mkdir(output_dir)
                
                
                % Example #2 for SBIEM:
                % An asperity breaks and triggers rupture in a neighboring asperity
                
                % get default parameters, with a few modified values:
                % play with a range of segment sizes (L) and resolutions (NX)
                % L allows control on the spurious periodicity effect
                % NX controls the quality of the solution
                % To explore efficiently the behavior of the asperities as a function
                % of model parameters it is convenient (shorter computing time) to
                % start with a small domain and low resolution.
                % It can be useful then to verify selected simulations by running
                % them again with 2*L and/or 2*NX
                
                % high resolution
                %pars = SBIEM([],'L',300e3,'NX',2048);
                
                % low resolution
                pars = SBIEM([],'L',400e3,'NX',1024*2*2);
                
                
                % high resolution, smaller domain:
                %pars = SBIEM([],'L',100e3,'NX',1024);
                
                pars.TMAX = 1.4 * pars.L/pars.CS;
                pars.OX_IX = [1:pars.NX];
                pars.OX_IDT = 4;
                pars.OT_IX = pars.NX/2-[0 40 80 120];
                
                
                
                
                
                % nucleation patch at maxima
                R_nuc=30e3; % range for nucleation
                
                max_TAU0=-1;
                TAU0=ones(size(pars.X))*100;
                XNUC=1;
                
                while abs(pars.X(XNUC))>R_nuc %(max_TAU0-pars.FRIC.MUd)<=0.8*(max(TAU0)-pars.FRIC.MUd) %
                    Phs=2*pi*rand(1,1024*2);
                    
                    % static and dynamic strength:
                    taus = pars.SIG0*pars.FRIC.MUs;
                    taud = pars.SIG0*pars.FRIC.MUd;
                    
                    % produce prestress distribution
                    
                    rough_model.name='Power';
                    rough_model.corrL=iCorrL*1e3;
                    rough_model.gamma=iHurst;
                    [~,TAU0]=rough_surfaceV2([-1 1]*pars.L/2,pars.NX,0.9.*(taus-taud),rough_model,Phs);
                    %TAU0=TAU0-min(TAU0)/4;
                    TAU0_1=TAU0;
                    
                    TAU0(TAU0<-0.0*(taus-taud))=-0.0*(taus-taud);
                    
                    temp_taper=tukeywin(pars.NX,0.9);
                    TAU0=TAU0.*temp_taper';
                    
                    
                    rough_model_BG.name='VonKarman';
                    rough_model_BG.corrL=1e3;
                    rough_model_BG.gamma=0.01;
                    % random background
                    [~,TAU0BG]=rough_surfaceV2([-1 1]*pars.L/2,pars.NX,0.05*(taus-taud),rough_model_BG,Phs_BG);
                    
                    %%
                    LBAR = max(pars.X)*0.5;%160e3 ;	% distance from center to barrier
                    TAU0=TAU0'+TAU0BG';
                    TAU0=TAU0/max(TAU0(abs(pars.X)<LBAR))*0.9.*(taus-taud);
                    
                    TAU0=TAU0+taud;
                    
                    % set low initial stress far out to stop rupture
                    SBAR = taud -0.5*(taus-taud) ; 	% stress in the barrier
                    
                    SBAR=0;
                    
                    TAU0(abs(pars.X)>LBAR)=taud;
                    
                    
                    %tau0( abs(X)>LBAR ) = SBAR;
                    taper_stress=tukeywin(length(pars.X),length(find(abs(pars.X)>LBAR))/length(pars.X));
                    TAU0=(TAU0-SBAR).*taper_stress+SBAR;
                                  
                    [max_TAU0,XNUC]=max(TAU0);
                    
                    
                end
                
                
                
                
                
                pars.TAU0=TAU0;
                
                if Index_test == 1
                    % plot the initial stress
                    figure
                    plot(pars.X,pars.TAU0/1e6,'LineWidth',1)
                    hold on
                    plot([pars.X(1) pars.X(end)],[taud taud]/1e6,'-m')
                    plot([pars.X(1) pars.X(end)],[taus taus]/1e6,'-b')
                    
                    pause
                    
                    continue
                end
                
                
                
                %% Nucleation length
                pars.FRIC.Dc=Dc_list(i_Dc);
                
                
                W = pars.SIG0.*(pars.FRIC.MUs-pars.FRIC.MUd)./ pars.FRIC.Dc; % slip weakening rate
                %if length(W)>1, W=W(NX/2); end % value at center
                W = max(W);
                LC = 2*0.57888694*pars.MU/W;
                
                disp(['Nucleation length: ' num2str(LC) ' m'])
                disp(['LC / DX = ' num2str(LC/pars.DX) ])
                
                %%  Nucleation around the middle at the local maxima
                
                LNUC = 1.5*LC;
                SNUC = pars.SIG0*(pars.FRIC.MUd +(pars.FRIC.MUs-pars.FRIC.MUd)*1.005);
                
                pars.TAU0(abs(pars.X-pars.X(XNUC))<=LNUC/2 ) = SNUC; % nucleation
                
                
                pars.CorrL=iCorrL;
                pars.Hurst=iHurst;
                
                
                temp_stress=(pars.TAU0-taud);
                temp_stress(temp_stress<0)=0;
                pars.avg_stressdrop=mean(temp_stress)/1e6;
                
                
                output_name=['Phs_' int2str(i_Phs) '_Dc_' num2str(Dc_list(i_Dc)) '_PW_Hurst_' num2str(iHurst)];
                
                disp(output_name)
                
                
                % do the computations
                [pars,xdat,tdat] = SBIEM(pars);
                vr = speed(xdat.X,xdat.RuptureTime);
                [T_extend,X_extend]=get_rupture_extends(xdat);
                
                
                %% plot the results
                close all
                disp('Plotting ...')
                
                FIGURE_name=[output_dir '/' output_name '.png'];
                figure(1)
                clf
                plot_simulations(xdat,tdat,pars,1,FIGURE_name);
                
                
                
                
                
                parsave([output_dir '/' output_name '.mat'],pars,xdat,tdat,taud,taus,vr)
                toc
                
                
                
                
                
            end
        end
    end
    
end


