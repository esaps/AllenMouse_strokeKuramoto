function[theta,  thds, i]=alloc_cnctmPlosCB(incoh, tau_m, varargin)
% allocate the history of phases,
% incoh is the proportion of unsynchronized oscillators
% if incoh >1, then oscilllators are arranged in antiphase
global mu h N Nt Nds

theta=zeros(N, tau_m);
thds=NaN(N, Nt);
%---
if incoh<=1
    theta(1:floor(incoh*N/2),        1) = linspace(pi/floor(incoh*N/2),     pi,floor(incoh*N/2));
    theta(end-floor(incoh*N/2)+1:end,1) = linspace(pi/floor(incoh*N/2)+pi,2*pi,floor(incoh*N/2));
else
    incoh=incoh-floor(incoh);
    theta(1:floor(N/2),1) = 1; % make the first half at pi
    theta(1:floor(incoh*N/2),        1) = linspace(pi/floor(incoh*N/2),     pi,floor(incoh*N/2));
    theta(floor(N/2)+1:end,1) = theta(floor(N/2),1)+pi;
end
thds(:,1)=theta(:,1);
%---------------------------------------------------------
% allocate the buffer from the delays
tt=0:h:(tau_m-h)*h;
for i=1:tau_m
    theta(:,i)=theta(:,1)+mean(mu)*tt(i);
    chck=ceil(i/Nds);
    if chck==i/Nds,
        thds(:,chck+1)=angle(mean(exp(1j*theta(:,(chck-1)*Nds+1:i)),2));
    end;
end
end