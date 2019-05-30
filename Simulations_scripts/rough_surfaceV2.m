function [x,surface_S1]=rough_surfaceV2(X,Nx,sigma,rough_model,fix_phase)
%% To produce parameterized rough surface with random phase perturbation
% [x,surface_S1]=rough_surface(X,Nx,sigma,rough_model,fix_phase)
% X: space range [X1 X2]
% Nx: number of points
% rough_model is a structure containing the information of rough model:
%
% for rough_model.name is 'Brune','VonKarman','DoubleCorrL', need to define
% rough_model.corrL (length scale of perturbation)
% rough_model.gamma (decay of high frequency components)
% fix_phase can be input as predefined for these cases
%
% rough_model.name can also be 'Predefined'
% in this situation, sigma can be any value and the fix_phase is actually
% the pre-defined distribution



if nargin==5
    Phs=fix_phase;
elseif nargin==4
    Phs=2*pi*rand(1,Nx/2); % random phase   
end




x=linspace(X(1),X(2),Nx); % m
dx=x(2)-x(1);
sk=1/dx/2;
k=(0:1:(Nx/2-1))/(Nx/2)*sk; % wavenumber
%rough_model='Brune';

gamma=rough_model.gamma;  % 
corr_length=rough_model.corrL;

switch rough_model.name
    
    case 'Predefined'
        surface_S1=fix_phase;
        return
        
    case 'Brune'
        %gamma=2;
        PSD=1./(1+k.^2*corr_length^2).^gamma;

    case 'VonKarman'
        %gamma=0.15;
        PSD=4*pi*corr_length^2*gamma./(1+k.^2*corr_length^2).^(1+gamma);
        PSD=PSD/max(PSD);
        
    case 'DoubleCorrL'
        %gamma=0.4;
        PSD=1./((1+k.^2*corr_length(1)^2).*(1+k.^2*corr_length(2)^2)).^(gamma); %1
        
    case 'Power'
        %gamma=0.4;
        PSD=10^(-4)*(2*pi)^3./k.^(gamma);
                          
end
%         figure(1)
%         loglog(k,PSD,'-o')
%         hold on


SP=sqrt(PSD).*exp(1j.*Phs); % spectrum


SP(end+1:Nx)=0;
SP(1)=SP(2);

surface_S=ifft(SP*Nx/2,'symmetric');
surface_S1=surface_S/max(abs(surface_S))*sigma;


% figure(5)
% subplot(211)
% hold on
% plot(x,surface_S1)

