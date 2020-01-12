close all
clear all
clc

% G5
% K4.3
% f4 [3-5] Hz
% G5K4.3f4v3
% diffKM86AMf2t1080K4.3D1h0.05strk0.9kp-0.25v3G5
% diffKM86AMf2t1080K4.3D1h0.05strk0.2kp-1v3Gp
% corr17regs4f_nsig_24w1G5
% corr17regs4f_nsig_24w4G5
% phcoh17cntr_rehab24f4 WTavg_control

% load(fullfile ('mouse_cnctm','CTXpl_iso_surface.mat'));
%
% for loops of different impacts of the stroke (Kstrk) and rebounds (keepK)
% in order to save time, the stroke is performed only once (one simulation
% of the time before stroke) phases before the stroke are saved (thds0 and
% theta0), as well as the integration step i0 when this occured, and later
% they are called so that simulation just continues for different levels of
% Kstrk and keepK
%%

Ares=25;
Aminvox=1;
Aminvol='1';
Afname=['res', int2str(Ares), 'minvox', int2str(Aminvox), 'minvol', Aminvol, '.h5'];
weights = h5read(Afname, '/weights');
lengths = h5read(Afname, '/tract_lengths');
reglabscell = h5read(Afname, '/region_labels'); % centres are the coordinates
centres = h5read(Afname, '/centres')'; % coordinates of the centres in mm
global N;
N=length(reglabscell);
reglabsn=NaN(N,1);
for i=1:N
    reglabsn(i)=length(reglabscell{i});
end
reglabs=repmat(' ',[N max(reglabsn)]);
for i=1:N
    reglabs(i, 1:reglabsn(i))=reglabscell{i};
end
ind = find(not(cellfun('isempty', strfind(reglabscell, 'olfactory'))));
Cxr=1:ind(1)-1;  Cxr=Cxr(:);
Cxl=Cxr+N/2; Cxl=Cxl(:);
Cx=[Cxr(:); Cxl(:)];
centres=centres/10; % unit is 100um, now it is 1mm
Cx_pos=centres(Cx, :);
Ncx=numel(Cx);
reglabscellCx=reglabscell(Cx);
%ds=0.2;
%load(['patchres25minvox1minvol1ds',num2str(ds),'zyx.mat']); % patches for the cortex
%load(fullfile ('mouse_cnctm','CTXpl_iso_surface.mat')); % surface of the cortex
load('mask24acr.mat');
load('regsNamesAcr.mat')
global mu fm K_cnctm h dt Nt Nds t D K0 tfin prntfgs hpar Nfmax incoh cnt0 Nh keepKs Kstrks
cxonly=1; % only the cortex
if cxonly==1
    N=Ncx;
end
keepKs=-1*[0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5];
Kstrks=[0:0.1:0.9];
if cxonly
    K0s=[0.5:0.1:8];
else
    K0s=[0.1:0.05:3];
end
Nr=24; % number of regions within the field of view
Ncond=4; % number of conditions, here 4 weeks compared to one week of control (always 5 measurements per week)
Igood={ ...
    [    6,    10,      13,14,      17   ]; ... % far from the upper left and completly in the box
    [  4,6,    10,      13,14,15,16,17   ]; ... % far from the upper left and completly in the box or have large part in the box
    [  4,6,    10,11,12,13,14,15,16,17   ]; ... % far from the upper left, minus regions which are mostly out of the box (strictest)
    [3,4,6,    10,11,12,13,14,15,16,17   ]; ... % far from the upper left, minus regions which are mostly out of the box (strict)
    [3,4,6,    10,11,12,13,14,15,16,17,18]; ... % far from the upper left, minus regions which are mostly out of the box (loose)
    [3,4,6,  9,10,11,12,13,14,15,16,17,18]; ... % far from the upper left, minus regions which are mostly out of the box (loosests)
    [3,4,6,8,9,10,11,12,13,14,15,16,17,18]; ... % almost every region apart MO
    [        9,10,11,12,13,14,15,16,17,18]; ... % far from the top
    [        9,10,      13,14,15,16,17,18]; ... % far from the upper left
    1:Nr};
printdiffG=1; % to print the difference of FC for different good regions; but also for detailed printing
fm=2;
D=1;
tfin=1080;
tadap=1;
if tadap, tfin0=tfin; end
incoh=0.8; % initial level of incoherence
T=40; % transient time to be discarded in the beggining and after the stroke
hpar=0.05;
Nfmax=5; % minimum 2; determines the step dt; dt=1/fmax=1/(Nfmax*fnyquist)
mu=2*pi*fm;
cnt0=500;
prntfgs=1;
Ngg=numel(Igood); % number of different choices of good nodes
%stroke=0.5; % 0.5 make stroke at half time; 0 no stroke, 1 stroke at the beggining
[difaS,difrS]              =deal(NaN(      numel(K0s),numel(keepKs), numel(Kstrks), Nr, Nr));

[cordifaS,cordifrS]        =deal(NaN(      numel(K0s),numel(keepKs), numel(Kstrks), Ncond, 15));
[cordifaGS,cordifrGS]      =deal(NaN(      numel(K0s),numel(keepKs), numel(Kstrks), Ncond, 15, Ngg));
txtname=['KM' int2str(N) 'AMf' num2str(fm) 't' num2str(tfin) 'a' num2str(tadap) 'D' num2str(D)...
    'K' num2str(K0s(1)) '_' num2str(K0s(end)) '_' int2str(numel(K0s)) ...
    'strk' num2str(Kstrks(1)) '_' num2str(Kstrks(end)) '_' int2str(numel(Kstrks)) ...
    'kp' num2str(keepKs(1)) '_' num2str(keepKs(end)) '_' num2str(numel(keepKs)) 'v2'];
disp(['f=', num2str(fm), ', D=', num2str(D), ', tfin=', num2str(tfin), ', N=', num2str(N)])
cnttmp=0;
for K0= K0s
    if tadap
        if K0>max(K0s)/2
            tfin=tfin0;
        else
            tfin=floor(tfin0 + 4*tfin0*(1 - (K0-min(K0s))/(max(K0s)/2-min(K0s)))); % 5*tfin0 for min(K0s) and tfin0 for max(K0s)/2
        end
    end
    IK0=find(K0==K0s);
    for keepK= keepKs
        IkeepK=find(keepK==keepKs);
        tic;
        for Kstrk=Kstrks
            IKstrk=find(Kstrk==Kstrks);
            disp(['K0=', num2str(K0), ', keepK=', num2str(keepK), ', Kstrk=', num2str(Kstrk)])
            rng('default'); rng(42);
            %%
            K_cnctm=K0*weights; % K_cnctm confirmed to the original at every loop
            if cxonly==1
                K_cnctm=K_cnctm(Cx, Cx);
            end
            %
            h=hpar/(max([max(K_cnctm(:)), mu*5*hpar, 5*D])); % at least 5 points per cycle. , max(sum(K_cnctm))
            dt=floor(pi/Nfmax/mu/h)*h; % fNyqust=2*fmax=Nfmax*mu/pi; dt=pi/Nfmax/mu~3/Nfmax/mu; for N=3.14 times larger frequency than f, 1 is obtained
            if dt<h, dt=h; end
            disp(['h=', num2str(h), '; dt=', num2str(dt), '; Kmax=', num2str(max(K_cnctm(:)), '%2.2f'), ...
                '; Kmean=', num2str(mean(K_cnctm(find(triu(ones(N),1)))), '%3.3f')])
            t=0:dt:tfin+2*dt; % so that it is ensured that tmax>=tfin; otherwise, if tfin+dt, with floor the last rounding (downsampling) might finish before tfin
            Nt=length(t);
            Nh=floor((tfin+2*dt)/h); %
            Nds=dt/h; % downsampling points, % there is a problem if h is not a multiple of dt
            if fix(Nds)~=Nds && abs(round(Nds)-Nds>1e-12),
                error('h should be a multiple of dt.');
            else
                Nds=round(Nds);
            end
            %%
            Tt=round(T/dt);
            eta=sqrt(2*D*dt)*rand(Ncx, Nt-Tt+1);
            if IkeepK==1 && IKstrk==1
                [thds, zs, ~, ~, tstrchck, theta0, i0]=KMcnctmAM(tfin/2, keepK, Kstrk);
            else
                [thds, zs, ~, ~, tstrchck]=KMcnctmAM(tfin/2, keepK, Kstrk, thds(:,1:tstrchck), theta0, i0);
            end
            %
            A=tril(ones(N), -1);
            [row, col]=find(A>0); clear A;
            Z=deal(NaN(N,N,2));
            for indregs1=1:N
                for indregs2=indregs1+1:N
                    Z(indregs1, indregs2, 1)=mean(exp(1j* ( thds(indregs2,Tt:floor(Nt/2)      )-thds(indregs1,Tt:floor(Nt/2)        ) ) ));
                    Z(indregs1, indregs2, 2)=mean(exp(1j* ( thds(indregs2,Tt+floor(Nt/2)+1:end)-thds(indregs1,Tt+floor(Nt/2)+1:end) ) ));
                    Z(indregs2, indregs1, :)=-Z(indregs1, indregs2, :);
                end
            end
            R=abs(Z);
            phcoh=angle(Z);
            fname=['KM' int2str(N) 'AMf' num2str(fm) 't' num2str(tfin) 'K' num2str(K0) 'D' num2str(D)...
                'h' num2str(hpar)  'strk' num2str(Kstrk) 'kp' num2str(keepK) 'v2']; % v2 after 20.12.17
            if (IkeepK==1 || IkeepK==numel(keepKs)) && (IKstrk==1 || IKstrk==numel(Kstrks))
                figure('units', 'centimeters', 'Position',[3 2 12*size(R,3) 11],'color','w','renderer','opengl')
                for i=1:size(R,3)
                    subplot(1,size(R,3),i)
                    Rt=squeeze(R(:,:,i));
                    him=imagesc(Rt(mask24acr, mask24acr));
                    set(him,'alphadata',~isnan(Rt(mask24acr, mask24acr)))
                    axis equal
                    set(gca, 'XTick', 1:Nr, 'XTickLabel', regexprep(acr270(mask24acr),{'_'},{''}), 'YTick', 1:Nr, 'YTickLabel', regexprep(acr270(mask24acr),{'_'},{''}), 'XTickLabelRotation', 45)
                    xlim([0.5, Nr+0.5]); ylim([0.5, Nr+0.5]);
                    clbar=colorbar;
                    caxis([min(min(min(R(mask24acr, mask24acr,:)))), max(max(max(R(mask24acr, mask24acr,:))))]);
                    clbar.Label.String = 'PLV';
                    if i==1
                        title(['Corr (sim) ', fname]);
                    else
                        title(['After stroke']);
                    end
                end
                try
                    saveas(gcf, fullfile('images', strcat(fname, '.png')), 'png');
                catch err
                end
            end
            %%
            difa= squeeze(R(mask24acr,mask24acr,2))-squeeze(R(mask24acr,mask24acr,1));
            difr=(squeeze(R(mask24acr,mask24acr,2))-squeeze(R(mask24acr,mask24acr,1))) ./ squeeze(R(mask24acr,mask24acr,1));
            difaS(IK0,IkeepK, IKstrk, :, :)= difa;
            difrS(IK0,IkeepK, IKstrk, :, :)= difr;
            cnttmp=cnttmp+1;
            disp([num2str(cnttmp/(numel(Kstrks)*numel(keepKs)*numel(K0s))*100), '%'])
            %                 end
        end
        pause(0.1); close all; pause(0.1);
        disp(['calculation time for one loop of Kstrks: ', num2str(toc/3600), ' hours, or ' , num2str(toc/60), ' minutes; ', int2str(numel(keepKs)*numel(K0s)), 'loops in total'])
    end
    disp(['calculation time for one loop of Kstrks and keepKs: ', num2str(toc/3600), ' hours, or ' , num2str(toc/60), ' minutes; ', int2str(numel(K0s)), 'loops in total'])
end
save(['mod',int2str(Nr), txtname, '.mat'], 'difaS', 'difrS', 'keepKs', 'K0s', 'Kstrks', '-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%