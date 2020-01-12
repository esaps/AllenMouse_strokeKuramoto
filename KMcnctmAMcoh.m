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
% to be added Igood = [(3),4,(6), 9, 10, 13, 14, 15, 16, 17]
%
% version 1 (no version) calculates the difference for each band as
% averaged over all points in that band, so it is skewed towards lower
% frequencies
% v2 calculates the difference as an integral under the curve for the
% frequency band, divided by the length of the band
%
% for loops of different impacts of the stroke (Kstrk) and rebounds (keepK)
% in order to save time, the stroke is performed only once (one simulation
% of the time before stroke) phases before the stroke are saved (thds0 and
% theta0), as well as the integration step i0 when this occured, and later 
% they are called so that simulation just continues for different levels of
% Kstrk and keepK


%%
printstat_fitt=0;
calcPrintstat=1;
savestat=0;
calcsim=0;
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
load('mask23-24acr.mat');
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
Nr=24;
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
%D=D*fm;
Dobserv=5;
tfin=1080;
tadap=1;
if tadap, tfin0=tfin; end
incoh=0.8;
T=40;
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
if calcsim
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
                    saveas(gcf, fullfile('images', strcat(fname, '.png')), 'png');
                end
                %%
                %                 if stroke==0.5
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
    pause(0.1); close all; pause(0.1);
    save(['mod',int2str(Nr), txtname, '.mat'], 'difaS', 'difrS', 'keepKs', 'K0s', 'Kstrks', '-v7.3');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print and calculate the statistics
if calcPrintstat
    Isupdiag=find(triu(ones(Nr),  1)>0);
    Isubdiag=find(tril(ones(Nr), -1)>0);
    [row, col]=find(tril(ones(Nr), -1)>0);
    Isubinverse= sub2ind([Nr,Nr], col, row);
    load(['mod',int2str(Nr), txtname, '.mat'])
    cnttmp=0;
    for K0= 4.3 %K0s
        IK0=find(K0==K0s);
        for keepK= -0.25% keepKs
            tic;
            IkeepK=find(keepK==keepKs);
            for Kstrk=0.9 %Kstrks
                IKstrk=find(Kstrk==Kstrks);
                % simulated data
                difa=squeeze(difaS(IK0,IkeepK, IKstrk, :, :));
                difr=squeeze(difrS(IK0,IkeepK, IKstrk, :, :));
                %
                eval(['mask=mask', int2str(Nr),'acr;'])
                w=weights(mask24acr, mask24acr);
                w(logical(eye(size(w))))=0;
                w=w/max(w(:)); % w is called in getprintstatFC; faster to be called outside than inside to be used load
                statsSim=NaN(3,Nr,Nr);
                statsSim(1,:,:)=w;
                statsSim(2,:,:)=difa;
                statsSim(3,:,:)=difr;
                %
                fname=['KM' int2str(N) 'AMf' num2str(fm) 't' num2str(tfin) 'K' num2str(K0) 'D' num2str(D)...
                    'h' num2str(hpar)  'strk' num2str(Kstrk) 'kp' num2str(keepK) 'v3'];
                %  save(fullfile('data_sim', [fname, 'dif.mat']), 'difa', 'difr', '-v7.3');
                for gg= 5 %1:numel(Igood) % print separately for each Igood
                    if gg==Ngg || printdiffG
                        figure('units', 'centimeters', 'Position',[3 2 24 8],'color','w','renderer','opengl')
                        for i=2:3
                            statstmp=squeeze(statsSim(i,:,:));
                            statstmp=statstmp(Igood{gg}, Igood{gg});
                            subplot(1,2,i-1); hold on;
                            set(gca, 'FontSize', 12)
                            him=imagesc(statstmp);
                            set(him,'alphadata',~isnan(statstmp))
                            axis square;
                            xlim([0.5, numel(Igood{gg})+0.5]); ylim([0.5, numel(Igood{gg})+0.5]);
                            if Nr==23 || Nr==24
                                set(gca, 'XTick', 1:numel(Igood{gg}), 'XTickLabel', regexprep(acr270(mask24acr(Igood{gg})),{'_'},{''}), ...
                                    'YTick', 1:numel(Igood{gg}), 'YTickLabel', regexprep(acr270(mask24acr(Igood{gg})),{'_'},{''}), 'XTickLabelRotation', 60, 'fontsize', 9)
                            else
                                set(gca, 'XTick', 1:numel(Igood{gg}), 'XTickLabel', maskacr(Igood{gg}), 'YTick', 1:numel(Igood{gg}), 'YTickLabel', maskacr(Igood{gg}), 'XTickLabelRotation', 45)
                            end
                            colorbar;
                            minclr=min(min(statstmp));
                            maxclr=max(max(statstmp));
                            caxis([minclr, maxclr])
                            %
                            if i==1
                                colormap(gca, 'parula')
                            else
                                % Create colormap that is green for negative, red for positive,
                                % and a chunk in the middle that is black.
                                Nclrg=256; Nclrr=Nclrg; % number of shades for each color
                                if abs(minclr)>abs(maxclr)
                                    Nclrg=round(Nclrg * abs(maxclr)/abs(minclr));
                                else
                                    Nclrr=round(Nclrr * abs(minclr)/abs(maxclr));
                                end
                                % also made colors symmetric, i.e. if in one dirrections values
                                % go "higher", make that color to become brighter, i.e.
                                % brightness to be absolute.
                                % also colors don;t start from black but from dark shades in
                                % order to get better differentiation
                                greenColorMap = [zeros(1, Nclrr), 0.15+linspace(0, 0.85*Nclrg/256, Nclrg)];
                                redColorMap = [linspace(0.85*Nclrr/256, 0, Nclrr)+0.15, zeros(1, Nclrg)];
                                colorMap = [redColorMap; greenColorMap; zeros(1, Nclrg+Nclrr)]';
                                % Apply the colormap.
                                colormap(gca, colorMap);
                            end
                            %
                            if i==1
                                if Nr==23 || Nr==24
                                    title(['SC weights and h=1 regions for ', int2str(Nr), ' regs'], 'fontsize', 10);
                                else
                                    title([ph(i), '-value for correlation for ', int2str(Nr), ' regs']);
                                end
                            elseif i==2
                                title(['abs. diff. N=', num2str(N), ', K0=', num2str(K0)], 'fontsize', 10);
                            elseif i==3
                                title(['rel. diff. keepK=', num2str(keepK), ', Kstrk=', num2str(Kstrk)], 'fontsize', 10);
                            end
                        end
                        if gg<Ngg
                            try
                                % /Volumes/Data/Spase's Data/mouse_connectivity/simKM
                                % /Volumes/Data/Spase's Data/Dropbox/matlab/ins/cnctms/mouse_cnctm/florence data
                                saveas(gcf, fullfile('images', 'diff', ['diff', fname, 'G', int2str(mod(gg, Ngg+1)), '.png']), 'png');
                            catch err
                                warning('error in finding the folder')
                            end
                        else
                            if (IkeepK==1 || IkeepK==numel(keepKs)) && (IKstrk==1 || IKstrk==numel(Kstrks))
                                saveas(gcf, fullfile('images', 'simKM', 'diff', ['diff', fname, '.png']), 'png');
                            else
                                try
                                    saveas(gcf, fullfile('images', 'diff', 'allregs',...
                                        ['diff', fname, '.png']), 'png');
                                catch err
                                    warning('error in finding the folder')
                                end
                            end
                        end
                    end
                end
                %  pause(0.1); close all; pause(0.1);
                %%
                % calculate statistics and print empirical data statistics
                % and FC
                % load(['pcorr',int2str(Nr),'regs4w0.mat']);
                tlimit=0;
                Nfb=7;
                load(['pcorr',int2str(Nr),'regs4w', int2str(round(tlimit)), 'Nb', int2str(Nfb), '.mat'])
                load(['phcoh17_',int2str(Nr),'regs', int2str(Nfb), 'fv2.mat']) % v2 corresponds to trapz integration for the frequency bands
                [cordifa   ,cordifr   ] = deal(NaN(Ncond, (Nfb + (Nfb>1))));
                [cordifaG  ,cordifrG  ] = deal(NaN(Ncond, (Nfb + (Nfb>1)), Ngg));
                [corPCdifa ,corPCdifr ] = deal(NaN(Ncond,  Nfb));
                [corPCdifaG,corPCdifrG] = deal(NaN(Ncond,  Nfb,            Ngg));
                IpermT = 1:Nr*(Nr-1)/2;
                stats=NaN(Nfb + (Nfb>1),Ncond,5,Nr,Nr);
                for ff=1:(Nfb + (Nfb>1))
                    % the last entry is for raw signals
                    if ss==1 % needs to be done onle once
                        [stats(ff,:,:,:,:)]=getprintstatFC(squeeze(pcorr17rehab_w1_4(:,:,ff,:)), squeeze(pcorr17control(:,:,ff)), ...
                            Nr, 'corr17regs', 0, 0, w); % gets statistics for each frequency band
                    end
                    % works for surrugoates and for actual data
                    
                    %
                    [difaT,difrT]=deal(NaN(Nr,Nr));
                    difaT(Isupdiag)=difa(Isupdiag(IpermT));
                    difaT(Isubdiag)=difaT(Isubinverse);
                    difrT(Isupdiag)=difr(Isupdiag(IpermT));
                    difrT(Isubdiag)=difrT(Isubinverse);
                    %
                    for iw=1:Ncond
                        %  empirical data, averaged filtered bands (by removing WT components)
                        difaemp=(stats(ff,iw,3,:,:));
                        difremp=(stats(ff,iw,4,:,:));
                        cordifa(iw, ff)=corr(difaemp(~isnan(difaemp)), difaT(~isnan(difaemp))); % corr(difaemp(intersect(find(~isnan(difaemp)), Isupdiag)), difa(intersect(find(~isnan(difaemp)), Isupdiag))); in order to have only the upper triangle, however the results are the same!
                        cordifr(iw, ff)=corr(difremp(~isnan(difremp)), difrT(~isnan(difremp)));
                        %  empirical data, averaged WT within the bands
                        if ff<=Nfb
                            difaempf=squeeze(difaf(iw,ff,:,:));
                            difrempf=squeeze(difrf(iw,ff,:,:));
                            corPCdifa(iw, ff)=corr(difaempf(~isnan(difaempf)), difaT(~isnan(difaempf)));
                            corPCdifr(iw, ff)=corr(difrempf(~isnan(difrempf)), difrT(~isnan(difrempf)));
                        end
                    end
                    for gg=1:numel(Igood)
                        difaGood=difaT(Igood{gg},Igood{gg});
                        difrGood=difrT(Igood{gg},Igood{gg});
                        for iw=1:Ncond
                            difaemp=(stats(ff,iw,3,Igood{gg},Igood{gg}));
                            difremp=(stats(ff,iw,4,Igood{gg},Igood{gg}));
                            cordifaG(iw, ff, gg)=corr(difaemp(~isnan(difaemp)), difaGood(~isnan(difaemp)));
                            cordifrG(iw, ff, gg)=corr(difremp(~isnan(difremp)), difrGood(~isnan(difremp)));
                            %
                            if ff<=Nfb
                                difaempf=(difaf(iw,ff,Igood{gg},Igood{gg}));
                                difrempf=(difrf(iw,ff,Igood{gg},Igood{gg}));
                                corPCdifaG(iw, ff, gg)=corr(difaempf(~isnan(difaempf)), difaGood(~isnan(difaempf)));
                                corPCdifrG(iw, ff, gg)=corr(difrempf(~isnan(difrempf)), difrGood(~isnan(difrempf)));
                            end
                        end
                    end
                end
                
                cordifaS( IK0,IkeepK, IKstrk, :,:  )=[corPCdifa, cordifa];
                cordifrS( IK0,IkeepK, IKstrk, :,:  )=[corPCdifr, cordifr];
                cordifaGS(IK0,IkeepK, IKstrk, :,:,:)=[corPCdifaG,cordifaG];
                cordifrGS(IK0,IkeepK, IKstrk, :,:,:)=[corPCdifrG,cordifrG];
                
            end
            cnttmp=cnttmp+1;
            disp([num2str(cnttmp/(numel(keepKs)*numel(K0s))*100), '% of printing simulations'])
            disp(['calculation time for one loop of Kstrks and keepKs: ', num2str(toc), ' seconds, or ' , num2str(toc/60), ' minutes; ', int2str(numel(keepKs)*numel(K0s)), ' loops in total'])
            pause(0.1); close all; pause(0.1);
        end
    end
    if savestat
        save(['cd17r',int2str(Nr),'f', int2str(Nfb), 'g', int2str(numel(Igood)), txtname(1:end-4),                        'v2.mat' ], 'cordifaS',    'cordifrS',    'cordifrGS',    'cordifaGS',...
            'Igood', 'keepKs', 'Kstrks', 'K0s',  '-v7.3')
    end
end
%%
if printstat_fitt
    cnttmp=0;
    if ~calcsim
        f_bands=[0, 0.15, 1, 3, 5, 7.5, 9, 12.5];
        Nfb=numel(f_bands)-1;
        load(['cd17r',int2str(Nr),'f', int2str(Nfb), 'g', int2str(numel(Igood)), txtname, '.mat'])
    end
    % print cordif
    for K0=K0s
        IK0=find(K0==K0s);
        cortmp=[cordifaGS(IK0,:,:,:,:,:), cordifrGS(IK0,:,:,:,:,:)];
        cmin=min(cortmp(:));
        cmax=max(cortmp(:));
        for iw=1:Ncond
            for gg=1:Ngg
                figure('units', 'centimeters', 'Position',[3 2 35 22],'color','w','renderer','opengl')
                for ff=1:7
                    subplot(4, 8, ff+8*0);
                    imagesc(Kstrks, keepKs, squeeze(cordifaGS(IK0,:,:,iw, ff, gg))); axis square; % plots first dimension along y, and the second along x; cordifaGS(IK0,IkeepK, IKstrk, :,:,:); but still it is called as imagesc(x,y,C)
                    %         surf(Kstrks, keepKs, squeeze(cordifaS(:,:,iw, ff)));
                    %         view(2); axis square;
                    caxis([cmin, cmax]);  axis square;
                    set(gca, 'XTick', Kstrks(1:2:end), 'fontsize', 8)
                    ytkstp=1;
                    set(gca, 'ytick',-5:ytkstp*5/9:0, 'yticklabel', (keepKs(end:-1*ytkstp:1)))
                    if ff==1
                        title(['avg. WT  f= ' num2str(f_bands(ff)) ' ' num2str(f_bands(ff+1))])
                    elseif ff==4
                        title(['abs. diff.  f= ' num2str(f_bands(ff)) ' ' num2str(f_bands(ff+1))])
                    else
                        title(['f= ' num2str(f_bands(ff)) ' ' num2str(f_bands(ff+1))])
                    end
                    %
                    subplot(4, 8, ff+8);
                    imagesc(Kstrks, keepKs, squeeze(cordifaGS(IK0,:,:,iw, ff+Nfb, gg))); axis square;
                    %set(gca, 'XTick', 1:2:numel(Kstrks), 'XTickLabel', Kstrks(1:2:end), 'YTick', 1:2:numel(keepKs), 'YTickLabel', keepKs(1:2:end), 'fontsize', 8)
                    set(gca, 'XTick', Kstrks(1:2:end), 'fontsize', 8)
                    ytkstp=1;
                    set(gca, 'ytick',-5:ytkstp*5/9:0, 'yticklabel', (keepKs(end:-1*ytkstp:1)))
                    caxis([cmin, cmax]);
                    if ff==1
                        title('filt. WT')
                    elseif ff==2
                        title(['w#', num2str(iw)])
                    elseif ff==3
                        title(['good', num2str(gg)])
                    elseif ff==5
                        title([int2str(Igood{gg})])
                    elseif ff==7
                        ax = gca;
                        axpos1 = ax.Position;
                    end
                    %
                    subplot(4, 8, ff+8*2);
                    imagesc(Kstrks, keepKs, squeeze(cordifrGS(IK0,:,:,iw, ff, gg))); axis square;
                    set(gca, 'XTick', Kstrks(1:2:end), 'fontsize', 8)
                    ytkstp=1;
                    set(gca, 'ytick',-5:ytkstp*5/9:0, 'yticklabel', (keepKs(end:-1*ytkstp:1)))
                    %  set(gca, 'XTick', 1:2:numel(Kstrks), 'XTickLabel', Kstrks(1:2:end), 'YTick', 1:2:numel(keepKs), 'YTickLabel', keepKs(1:2:end), 'fontsize', 8)
                    caxis([cmin, cmax]);
                    if ff==1
                        title('avg. WT')
                    elseif ff==2
                        title(['K=', num2str(K0), ', N=', num2str(N)])
                    elseif ff==3
                        title(['f=', num2str(fm), ', D=', num2str(D)])
                    elseif ff==4
                        title('rel. diff.')
                    end
                    %
                    subplot(4, 8, ff+8*3);
                    imagesc(Kstrks, keepKs, squeeze(cordifrGS(IK0,:,:,iw, ff+Nfb, gg))); axis square;
                    set(gca, 'XTick', Kstrks(1:2:end), 'fontsize', 8)
                    ytkstp=1;
                    set(gca, 'ytick',-5:ytkstp*5/9:0, 'yticklabel', (keepKs(end:-1*ytkstp:1)))
                    %set(gca, 'XTick', 1:2:numel(Kstrks), 'XTickLabel', Kstrks(1:2:end), 'YTick', 1:2:numel(keepKs), 'YTickLabel', keepKs(1:2:end), 'fontsize', 8)
                    caxis([cmin, cmax]);
                    if ff==1
                        title('filt. WT')
                    end
                    %
                end
                subplot(4, 8, 8);
                imagesc(Kstrks, keepKs, squeeze(cordifaGS(IK0,:,:,iw, 8, gg)));
                set(gca, 'XTick', Kstrks(1:2:end), 'fontsize', 8)
                ytkstp=1;
                set(gca, 'ytick',-5:ytkstp*5/9:0, 'yticklabel', (keepKs(end:-1*ytkstp:1)))
                %  set(gca, 'XTick', 1:2:numel(Kstrks), 'XTickLabel', Kstrks(1:2:end), 'YTick', 1:2:numel(keepKs), 'YTickLabel', keepKs(1:2:end), 'fontsize', 8)
                caxis([cmin, cmax]);  axis square;
                title('raw')
                %
                subplot(4, 8, 8+8*2);
                imagesc(Kstrks, keepKs, squeeze(cordifrGS(IK0,:,:,iw, 8, gg)));
                set(gca, 'XTick', Kstrks(1:2:end), 'fontsize', 8)
                ytkstp=1;
                set(gca, 'ytick',-5:ytkstp*5/9:0, 'yticklabel', (keepKs(end:-1*ytkstp:1)))
                % set(gca, 'XTick', 1:2:numel(Kstrks), 'XTickLabel', Kstrks(1:2:end), 'YTick', 1:2:numel(keepKs), 'YTickLabel', keepKs(1:2:end), 'fontsize', 8)
                caxis([cmin, cmax]);  axis square;
                title('raw')
                % plot colorbar
                subplot(4, 8, 8+8);
                axis off;
                cbh=colorbar;
                cpos = cbh.Position;
                cl=cpos(4)-cpos(2);
                cpos(4) = axpos1(4);
                cpos(2) = axpos1(2);
                cpos(3) = 2*cpos(3);
                cbh.Position = cpos;
                caxis([cmin, cmax]);
                saveas(gcf, fullfile('images', 'ParamSpace', ['goodreg', '17r',int2str(Nr),'f', int2str(Nfb), 'KM' int2str(N), ...
                        'f' num2str(fm) 't' num2str(tfin) 'K' num2str(K0) 'D' num2str(D), 'w', int2str(iw), 'g', int2str(gg), 'v', int2str(2+ (keepKs(end)<0)),'.png']), 'png');
            end
            pause(0.1); close all; pause(0.1);
            cnttmp=cnttmp+1;
            disp([num2str(cnttmp/(numel(Ncond)*numel(K0s))*100), '% of printing simulations'])
        end
    end
end
%%
