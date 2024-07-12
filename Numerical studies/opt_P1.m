function [flag,w_1,net_f_P1,TSO_Goal,DSO_Goal] = opt_P1(M,N,lambda,tran_one,G,DD,agg_mat,eta,H_T,H_D,pp_k,S_max,s_max,all_mat,C_D,C_U)

flag=0;
w_0=(lambda*abs(tran_one*(G-DD))+ eta* agg_mat'*(max(0,abs(H_D*(pp_k))-s_max))*all_mat);
cvx_begin
    variable f(M,N);
    minimize(lambda*abs(tran_one*(G-DD+f*agg_mat))+ eta* agg_mat'*(max(0,abs(H_D*(pp_k-f'))-s_max))*all_mat);
    subject to
    f <= C_D;
    f >= C_U;
    H_T*(G-DD+f*agg_mat) <= S_max;
    H_T*(G-DD+f*agg_mat) >= -S_max;

cvx_end
net_f_P1=f;
TSO_Goal=abs(tran_one*(G-DD+f*agg_mat));
DSO_Goal=agg_mat'*(max(0,abs(H_D*(pp_k-f'))-s_max));
w_1=cvx_optbnd;
if w_1>w_0
    flag=1;
    net_f_P1=zeros(M,N);
    w_1=w_0;
    TSO_Goal=abs(tran_one*(G-DD));
    DSO_Goal=agg_mat'*(max(0,abs(H_D*(pp_k))-s_max));
end
end

