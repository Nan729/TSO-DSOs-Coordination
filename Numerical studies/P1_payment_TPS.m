function [T_to_D,f_get] = P1_payment_TPS(w_0,w_1,bar_pow,u_T,u_D,M,N,DSO_Goal)

alpha=bar_pow(1);
beta=bar_pow(2);
gamma=bar_pow(3);
f_get=gamma*(w_0-w_1);
u_T(1,2)=alpha*(w_0-w_1)+u_T(1,1);
u_D(1,2)=beta*(w_0-w_1)+u_D(1,1);
u_D(2,2)=beta*(w_0-w_1)+u_D(1,1);
mediate=u_D(:,2)+f_get+DSO_Goal';
T_to_D=mediate(1:2,:);
end
