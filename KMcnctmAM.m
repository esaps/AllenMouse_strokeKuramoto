%% 
function[thds, zs, zs1, zs2, tstrchck, theta0, i0]=KMcnctmAM(varargin)
% function for calculating phases and mean field of Kuramoto oscillators 
% for a connectome with white noise
% thds  are phases for all t (downsampled) and N nodes
% zs is the Kuramoto order parameter (zs1 and zs2 are KOP for each of the
% hemispheres, not important here)
% for loops of different impacts of the stroke (Kstrk) and rebounds (keepK)
% in order to save time, the stroke is performed only once (one simulation
% of the time before stroke) phases before the stroke are saved (thds0 and
% theta0), as well as the integration step i0 when this occured, and later 
% they are called so that simulation just continues for different levels of
% Kstrk and keepK

global mu K_cnctm h dt N Nt Nds t D tfin incoh cnt0 Nh par fm K0 hpar keepKs Kstrks
% allocation
stroke=0; tstrchck=NaN;
strkflag=1; % to allow one stroke if the other conditions are fullfilled
if nargin>0
    tstr = varargin{1}; % time of the stroke
    stroke=1; % if there is time of the stroke defined, it means there will be a stroke during the simulation, hence stroke=1;
    tstrchck=find(tstr<t, 1, 'first'); % position in the time-series of the time when the stroke occurs 
    keepK=1; % to keep same level of K, or to increase
    Kstrk=0.5; % how much to be removed
    if stroke && nargin>1
        keepK=varargin{2};
        if nargin>2
            Kstrk=varargin{3};
            if nargin>3
                thds0=varargin{4};
                theta0=varargin{5};
                i0=varargin{6};
            end
        end
    end
end
if nargin<6
    [theta, thds, i]=alloc_cnctmPlosCB(incoh, Nds+1);
else
    theta=theta0;
    i=i0; % to be careful not to cut the Kcnctm twice!
    thds(:,1:tstrchck)=thds0;
end
%%
varstrKM.mu=mu;
varstrKM.K_cnctm=K_cnctm;
varstrKM.h=h;
varstrKM.Nds=Nds;
varstrKM.t=t;
%  varstrKM.cnt=cnt0-mod(cnt0,dt);
varstrKM.cnt=round(cnt0/dt);
varstrKM.N=N;
varstrKM.N2=N/2;
while i<Nh
    eta=sqrt(2*h*D)*randn(N,1);
    chck=ceil((i+1)/Nds);   % counter for the downsampled phases until chck=(i+1)/Nds, all the thds and thds2 are returned NaN
    Ihealthy=setdiff(1:N,2);
    if stroke && chck==tstrchck && strkflag
        if keepK>=0
            K_cnctm(Ihealthy,Ihealthy)=K_cnctm(Ihealthy,Ihealthy) * (1 + keepK*(sum(K_cnctm(:,2))+sum(K_cnctm(2,:)))/sum(K_cnctm(:))); % sustain the same level of coupling (for keepK=1) or maybe compansate to some percent, e.g. keepK=0.5
        else
            % the links which had links to M1, will get links to eachother
            K2out=K_cnctm(2,Ihealthy)/sum(K_cnctm(2,Ihealthy)); % normalized outgoing links            
            for ii=Ihealthy
                if K_cnctm(ii,2)>0
                    K_cnctm(ii,Ihealthy)= K_cnctm(ii,Ihealthy) + K_cnctm(ii,2) * abs(keepK) * K2out; % incoming to 2
                    % each incoming link to 2, is switched equaly towards all the outgoing directions from 2. 
                end
            end
        end   
        K_cnctm(2,:)=Kstrk * K_cnctm(2,:);
        K_cnctm(:,2)=Kstrk * K_cnctm(:,2); % incoming to 2
        varstrKM.K_cnctm=K_cnctm;
        strkflag=0; % to avoid "stroking" more than once
        theta0=theta; i0=i; % to keep those for the next cycle
    end
    
    [theta, thds(:,chck+1), i] = KMcnctmHt0(theta, eta, i, chck, varstrKM);                             % returns i+1
end
%% ---------------
% still to be able to print in case of breaking the simulation
[~,I]=find(isnan(thds), 1, 'first');
if isempty(I),
    I=length(thds(1,:))+1;
else
    if I<length(t)-1
        warning('NaN detected|')
        disp(['I=', int2str(I), ',  t=', num2str(t(I))])
        pause;
    end
end
t(I:end)=[];
Nt=length(t);
tfin=floor(t(I-1)); % to keep it integer
thds(:,I:end,:)=[];
thds=mod(rem(thds+pi, pi*2)+pi*2, pi*2)-pi;
%
par.fnamedir='AMstroke';
fname=['AM' int2str(N) 'f' num2str(fm) 't' num2str(tfin) 'K' num2str(K0) 'D' num2str(D)...
    'h' num2str(hpar)  'strk' num2str(Kstrk) 'kp' num2str(keepK)];
% try
%     if (find(keepK==keepKs)==1 || find(keepK==keepKs)==numel(keepKs)) && (find(Kstrk==Kstrks)==1 || find(Kstrk==Kstrks)==numel(Kstrks))
%         [zs, zs1, zs2]=printFTevol(thds, fname, 1);
%     else
%         [zs, zs1, zs2]=printFTevol(thds, fname, 0);
%     end
% catch err
%     [zs, zs1, zs2]=printFTevol(thds, fname, 0);
% end
zs = squeeze(mean(exp(1j*thds),1));
[zs1, zs2]=deal(NaN);

end

