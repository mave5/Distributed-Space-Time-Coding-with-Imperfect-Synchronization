% this function coverts MPSK symbols to binary 0/1 
% gray mapping is used here 
% M=2 0->-1 , 1-> +1
% M=4 00->1,01->j,10->-j,11->-1

% M=8, 8-PSK
% 000 -> -1-j
% 001 -> -1
% 010 -> j
% 011 -> -1+j
% 100 -> -j
% 101 -> 1-j
% 110 -> 1+j
% 111 -> 1

function bhat=mpsk2bin(s,M)
% s : Mpsk symbols, the symbols can be noisy or free of noise
% M : M-psk constelation size
% bhat: binary 0/1 data as a vector
 
N=length(s);
M8_map=[1,-1i,1i,-1,(1+1i)/sqrt(2),(-1+1i)/sqrt(2),(-1-1i)/sqrt(2),(1-1i)/sqrt(2)];
bhat=zeros(log2(M)*N,1);

if M==2 
    bhat(:,1)=(sign(real(s))+1)/2;
elseif M==4    
% rotate pi/4, to make decision rules simpler
s= s*exp(1i*pi/4); 
%bhat=zeros(2*N,1);
bhat(1:2:end)=.5*(1-sign(imag(s)));
bhat(2:2:end)=.5*(1-sign(real(s)));
elseif M==8       
         for k=1:length(s)
             [~,ind]=min(abs(M8_map-s(k)));
            switch M8_map(ind)
                case 1
                    bhat(3*k-2:3*k)=[1;1;1];
                case -1i
                    bhat(3*k-2:3*k)=[1;0;0];
                case 1i
                    bhat(3*k-2:3*k)=[0;1;0];
                case -1
                    bhat(3*k-2:3*k)=[0;0;1];
                case (1+1i)/sqrt(2)
                    bhat(3*k-2:3*k)=[1;1;0];
                case (-1+1i)/sqrt(2)
                    bhat(3*k-2:3*k)=[0;1;1];
                case (-1-1i)/sqrt(2)
                    bhat(3*k-2:3*k)=[0;0;0];
                case (1-1i)/sqrt(2)
                    bhat(3*k-2:3*k)=[1;0;1];
            end    
        end
end