% " Differential Distributed Space-Time Coding with Imperfect Synchronization", IEEE GLOBECOM 2014, Austin, USA
% M. R. Avendi, 2014

close all;
clear all;
clc;
addpath functions;
%% parameters

% total number of symbols to be simulated
% if you want to make simulation fast, you can reduce it
% however to get a smooth plot, you need around 2E5 symbols
Nsyms=1E4;

% number of relays
R=2; 

tau1=0; % Relay 1 delay
tau2=0.5; % Relay 2 delay 2

% match filter parameters
a=.9;
alfa1=sinc(tau1).*cos(pi*a*tau1)./(1-4*a^2*tau1^2);
tau11=1-tau1;
beta11=sinc(tau11).*cos(pi*a*tau11)./(1-4*a^2*tau11^2);

% raised cosine matched filter
a=0.9;
alfa2=sinc(tau2).*cos(pi*a*tau2)./(1-4*a^2*tau2^2);
tau22=1-tau2;
beta2=sinc(tau22).*cos(pi*a*tau22)./(1-4*a^2*tau22^2);

M=2; % MPSK symbols
Ptot_dB=0:5:30;% totla power
N0=1;
Ptot=10.^(Ptot_dB/10)*N0;

% Alamouti type 
type=1;

% channel variances
sig_sr=1;

% power allocation and amplification factor in relays
P0=Ptot./2;
Pr=Ptot./4;
AF= Pr./(P0*sig_sr+N0);

%number of OFDM symbols
N=64;
Ncp=1;% cyclic prefix length

% make sure total number of symbols is a multiple of N
Ns=N*floor(Nsyms/N);

% channel correlation 
fdTs=1e-3;
% channel use 
ch_use=1;

%% mail loop 
for snr_ind=1:length(Ptot)

% counters for number of errors     
nerr1=0;
nerr2=0;

nbits=0;
nSim=0;

display('simulation in progress ....')
current_power_simulation=Ptot_dB(snr_ind)


% generate channels
h1=flat_cos(Ns,fdTs,ch_use); % Source to Relay channels
h2=flat_cos(Ns,fdTs,ch_use);

g1=flat_cos(Ns,fdTs,ch_use); % relay to destination channels
g2=flat_cos(Ns,fdTs,ch_use);

% if you want to make simulation faster, you can reduce Ns
while  nSim<Ns-N 

% symbol counter   
nSim=nSim+1;

%input bits generation
b1_in=bits(log2(M)*N);
b2_in=bits(log2(M)*N);

% MPSK symbol
v1_in=bin2mpsk(b1_in,M); 
v2_in=bin2mpsk(b2_in,M); 
temp=[v1_in,v2_in];
v_in=reshape(temp.',2*N,1);

% space-time encoding
V_in=stc_alamouti(v_in,type);

% differential encoder
if nSim==1
    s_km1=[ones(1,N);ones(1,N)]/sqrt(2);
    s_k=s_km1;
else
    s_k=diff_encoder_v(V_in,s_km1);
    s_km1=s_k;
end

% OFDM
S1=sqrt(N)*ifft(s_k(1,:));
S2=sqrt(N)*ifft((s_k(2,:)));

% AWGN noise CN(0,N0)
n11=cxn(N,N0);%noise at relay 1, TS1
n12=cxn(N,N0);%noise at relay 1, TS2
n21=cxn(N,N0);%noise at relay 2, TS1
n22=cxn(N,N0);%noise at relay 2, TS2

% RX signals at relays
X11=sqrt(P0(snr_ind)*R)*h1(nSim)*S1+n11;
X12=sqrt(P0(snr_ind)*R)*h1(nSim)*S2+n12;
X21=sqrt(P0(snr_ind)*R)*h2(nSim)*S1+n21;
X22=sqrt(P0(snr_ind)*R)*h2(nSim)*S2+n22;

% circular time-reverse
X21_tr=[X21(1) fliplr(X21(2:length(X21)))];
X22_tr=[X22(1) fliplr(X22(2:length(X22)))];

% Add Cyclic Prefrex
CP11=X11(end-Ncp+1:end);
CP12=X12(end-Ncp+1:end);
CP21=X21_tr(end-Ncp+1:end);
CP22=X22_tr(end-Ncp+1:end);
X11_cp=sqrt(AF(snr_ind))*[CP11,X11];
X12_cp=sqrt(AF(snr_ind))*[CP12,X12];
X21_cp=sqrt(AF(snr_ind))*[CP21,X21_tr];
X22_cp=sqrt(AF(snr_ind))*[CP22,X22_tr];

% relay to destination channels 
g1t=[g1(nSim),zeros(1,N-1)];
g2t=[alfa2*g2(nSim),beta2*g2(nSim),zeros(1,N-2)];

% recevied signals at destination
Y1_ch=conv(g1t,X11_cp)-conv(g2t,conj(X22_cp)); % convolution
Y2_ch=conv(g1t,(X12_cp))+conv(g2t,conj(X21_cp));

% remove tails
Np=N+Ncp;
Y1_ch=Y1_ch(1:Np);
Y2_ch=Y2_ch(1:Np);

% add AWGN 
ly=length(Y1_ch);
w1=cxn(ly,N0);
w2=cxn(ly,N0);
Y1_cp=Y1_ch+w1;
Y2_cp=Y2_ch+w2;
 
%%  Decoding process at destiatnion
% remove CP
 Y1=Y1_cp(Ncp+1:end);
 Y2=Y2_cp(Ncp+1:end);
 
  % demodulation OFDM
 y1= sqrt(1/N)*fft(Y1,N);
 y2= sqrt(1/N)*fft(Y2,N);
 y_k=[y1;y2];
 
% non-coherent decoding
 if nSim==1
     y_km1=y_k;
 else
    [v1_hat,v2_hat]= dstc_decoder(y_k,y_km1,type);
    y_km1=y_k;
    
    % MPSK demodulation
    b1_hat=mpsk2bin(v1_hat,M);
    b2_hat=(mpsk2bin(v2_hat,M));
    
    % count errors
    nerr1=nerr1+sum(abs(b1_in-b1_hat));
    nerr2=nerr2+sum(abs(b2_in-b2_hat));
    nbits=log2(M)*N+nbits;
 end
 
end

% Bit Error Rate
BER1(snr_ind)=(nerr1)/(nbits);
BER2(snr_ind)=(nerr2)/(nbits);

end

%% PLOT 
Pb_th=1./(4*Ptot);
clr=['r-+'; 'k->'; 'y-*'; 'g-o'; 'r-<'];
figure
semilogy(Ptot_dB,BER1,clr(1,:),'LineWidth',2,'MarkerSize',8); 
hold on;
semilogy(Ptot_dB,BER2,clr(4,:),'LineWidth',2,'MarkerSize',8); 
grid on
xlabel('Total Power');
ylabel('BER');
legend(['D-OFDM DSTC, \tau=',num2str(tau2)]);

set(gca,'XTick',Ptot_dB(1):5:Ptot_dB(end),'FontSize',16,...
   'FontName','Times New Roman');
