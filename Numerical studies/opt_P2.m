function [flag,w_2,net_f_P2,TSO_Goal,DSO_Goal] = opt_P2(M,N,lambda,tran_one,G,DD,agg_mat,eta,H_T,H_D,pp_k,S_max,s_max,all_mat,C_D,C_U)

flag=0;
w_0=(lambda*abs(tran_one*(G-DD))+ eta* agg_mat'*(max(0,abs(H_D*(pp_k))-s_max))*all_mat);
cvx_begin
    variable f(M,N);
    minimize(lambda* abs(tran_one*(G-DD+f*agg_mat)));
    subject to
    f <= C_D;
    f >= C_U;
    H_T*(G-DD+f*agg_mat) <= S_max;
    H_T*(G-DD+f*agg_mat) >= -S_max;
    H_D*(pp_k-f') <= s_max;
    H_D*(pp_k-f') >= -s_max;
cvx_end
net_f_P2=f;
TSO_Goal=abs(tran_one*(G-DD+f*agg_mat));
DSO_Goal=zeros(1,M);
w_2=cvx_optbnd;
if w_2>w_0
    flag=1;
    net_f_P2=zeros(M,N);
    w_2=w_0;
    TSO_Goal=abs(tran_one*(G-DD));
    DSO_Goal=agg_mat'*(max(0,abs(H_D*(pp_k))-s_max));
end
end

