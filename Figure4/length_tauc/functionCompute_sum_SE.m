function [SE_sum] = functionCompute_sum_SE(SE_u,SE_d)


SE_sum=0.5*sum(SE_u(:)+0.5*SE_d(:));
