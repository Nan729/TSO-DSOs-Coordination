function [T_to_D,D_to_f] = P2_payment(w_0,w_2,bar_pow,u_T,u_D,M,N)

alpha=bar_pow(1);
beta=bar_pow(2);
gamma=bar_pow(3);
D_to_f=gamma*(w_0-w_2);
u_T(1,2)=alpha*(w_0-w_2)+u_T(1,1);
u_D(1,2)=beta*(w_0-w_2)+u_D(1,1);
u_D(2,2)=beta*(w_0-w_2)+u_D(1,1);
T_to_D=u_D(1:2,2)+D_to_f;
end

