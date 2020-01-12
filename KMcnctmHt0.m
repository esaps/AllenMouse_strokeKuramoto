%-------------------------------------------------------------------------%
%-                           Spase Petkoski                              -%
%-Theoretical Neuroscience Group--Institut de Neurosciences des Syst?mes -%
%-                        2019, AMU Marseille                            -%
%-------------------------------------------------------------------------%
% function for calculating phases and mean field at each step of the Heun
% integration of Kuramoto oscillators for a connectome with white noise
% theta are the phases up to history t-tau_m
% th are the downsampled phases, kept for all t and N nodes
% no time delays
function [theta, th, i] = KMcnctmHt0(theta, eta, i, chck, varstrKM)
%-------------
% global mu K_cnctm h N Nds t cnt
if isfield(varstrKM, 'omega')
    omega=varstrKM.omega;
else
    omega=varstrKM.mu;
end
K_cnctm=varstrKM.K_cnctm;
h=varstrKM.h;
Nds=varstrKM.Nds;
t=varstrKM.t;
cnt=varstrKM.cnt;
N=varstrKM.N;
%-------------------------------------------------------------------------
%------- connectome with no delays------- mean field approach ------------
%-------------------------------------------------------------------------
zi_t= (K_cnctm*exp(1j*theta(:,end)));
kn13 = omega - abs(zi_t).*sin(theta(:,end) - angle(zi_t));

zi_t= (K_cnctm*exp(1j*(theta(:,end) + h*kn13 + eta)));
kn23 = omega - abs(zi_t).*sin(theta(:,end) + h*kn13 + eta - angle(zi_t));

theta(:,1:end-1)=theta(:,2:end);
theta(:,    end)=theta(:,end) + h/2*(kn13+kn23) + eta;
%---
i=i+1;
%---
%chck=ceil(i/Nds); %after the condition below is true, move to next chck and all of them are NaN until the next integer i/Nds
if chck==i/Nds
    th=angle(mean(exp(1j*theta(:,end-Nds:end)),2));
    if ~mod(chck,cnt)
        disp(round(t(chck)));
    end
else
    th=NaN(N,1);
end
end