function [Phi,F] = mpcgain(A_aug,B_aug,C_aug,N_c,N_p);
% [m1,n1]=size(Cp);
% [n1,n_in]=size(Bp);
% A_e=eye(n1+m1,n1+m1);
% A_e(1:n1,1:n1)=Ap;
% A_e(n1+1:n1+m1,1:n1)=Cp*Ap;
% B_e=zeros(n1+m1,n_in);
% B_e(1:n1,:)=Bp;
% B_e(n1+1:n1+m1,:)=Cp*Bp+Dp;
% C_e=zeros(m1,n1+m1);
% C_e(:,n1+1:n1+m1)=eye(m1,m1);
% D_e = Dp;
% 
% n=n1+m1;
h(1,:)=C_aug;
F(1,:)=C_aug*A_aug;
for kk=2:N_p
h(kk,:)=h(kk-1,:)*A_aug;
F(kk,:)= F(kk-1,:)*A_aug;
end
v=h*B_aug;
Phi=zeros(N_p,N_c); %declare the dimension of Phi
Phi(:,1)=v; % first column of Phi
for i=2:N_c
Phi(:,i)=[zeros(i-1,1);v(1:N_p-i+1,1)]; %Toeplitz matrix
end
% BarRs=ones(N_p,1);
% Phi_Phi= Phi'*Phi;
% Phi_F= Phi'*F;
% Phi_R=Phi'*BarRs;