%% Compare jns.m to size & power from
% Ledoit & Wolf (2002). Some hypothesis tests for the covariance matrix
%   when the dimension is large compared to the sample size.
%   Annals of Statistics 30: 1081-1102

%% Table 1
clear all;
n = [4 8 16 32 64 128 256];
p = [4 8 16 32 64 128 256];
reps = 200;
tic;
for i = 1:numel(p)
   for j = 1:numel(n)
      for k = 1:reps
         x = randn(n(j),p(i));
         pval(k) = sphere.jsn(x,'test','wang');
      end
      prob(i,j) = mean(pval<=0.05);
   end
   toc
end

% Ledoit & Wolf (method is John's test)
pL = [...
0.01 0.03 0.04 0.05 0.05 0.05 0.05;...
0.03 0.04 0.04 0.05 0.05 0.05 0.05;...
0.04 0.05 0.05 0.05 0.05 0.05 0.05;...
0.05 0.05 0.05 0.05 0.05 0.05 0.05;...
0.05 0.05 0.05 0.05 0.05 0.05 0.05;...
0.05 0.05 0.05 0.05 0.05 0.05 0.05;...
0.05 0.05 0.05 0.05 0.05 0.05 0.05];

% method = 'john', reps = 2000; 25.11.14
% prob =
% 
%     0.0010    0.0230    0.0400    0.0500    0.0415    0.0505    0.0535
%     0.0220    0.0370    0.0435    0.0485    0.0520    0.0490    0.0635
%     0.0335    0.0415    0.0460    0.0515    0.0605    0.0460    0.0470
%     0.0445    0.0505    0.0490    0.0460    0.0505    0.0505    0.0460
%     0.0470    0.0515    0.0560    0.0535    0.0530    0.0480    0.0500
%     0.0565    0.0545    0.0505    0.0515    0.0520    0.0515    0.0605
%     0.0580    0.0510    0.0530    0.0530    0.0575    0.0550    0.0490

% method = 'nagao', reps = 2000; 25.11.14
% prob =
% 
%     0.0420    0.0400    0.0540    0.0485    0.0495    0.0455    0.0510
%     0.0720    0.0585    0.0540    0.0490    0.0475    0.0455    0.0550
%     0.0310    0.0510    0.0485    0.0410    0.0555    0.0530    0.0530
%     0.0280    0.0430    0.0430    0.0510    0.0475    0.0455    0.0505
%     0.0315    0.0475    0.0480    0.0485    0.0530    0.0450    0.0495
%     0.0370    0.0430    0.0410    0.0475    0.0565    0.0520    0.0535
%     0.0355    0.0415    0.0450    0.0495    0.0435    0.0525    0.0370

% method = 'wang', reps = 2000; 25.11.14
% prob =
% 
%     0.0270    0.0285    0.0445    0.0495    0.0485    0.0515    0.0620
%     0.0505    0.0595    0.0380    0.0495    0.0485    0.0575    0.0535
%     0.0495    0.0750    0.0635    0.0490    0.0695    0.0555    0.0570
%     0.0500    0.0780    0.0965    0.0780    0.0575    0.0585    0.0550
%     0.0475    0.0900    0.1030    0.1090    0.0725    0.0485    0.0480
%     0.0525    0.0860    0.1030    0.1125    0.1215    0.0845    0.0455
%     0.0585    0.0795    0.1070    0.1135    0.1250    0.1140    0.0885

%% Table 2
clear all;
n = [4 8 16 32 64 128 256];
p = [4 8 16 32 64 128 256];
reps = 200;
tic;
for i = 1:numel(p)
   for j = 1:numel(n)
      for k = 1:reps
         sigma = [0.5*ones(round(p(i)/2),1);ones(p(i)-round(p(i)/2),1)]';
         x = (diag(sqrt(sigma))*randn(n(j),p(i))')';
         pval(k) = sphere.jsn(x,'test','w');
      end
      prob(i,j) = mean(pval<=0.05);
   end
   toc
end

% Ledoit & Wolf (method is John's test)
pL = [...
0.02 0.06 0.15 0.37 0.76 0.98 1;...
0.05 0.09 0.18 0.42 0.85 1.00 1;...
0.06 0.11 0.20 0.48 0.90 1.00 1;...
0.08 0.13 0.22 0.50 0.93 1.00 1;...
0.09 0.13 0.24 0.52 0.95 1.00 1;...
0.09 0.14 0.23 0.53 0.95 1.00 1;...
0.09 0.14 0.24 0.54 0.96 1.00 1];

% method = 'john', reps = 2000; 25.11.14
% prob =
% 
%     0.0010    0.0505    0.1310    0.3425    0.7420    0.9880    1.0000
%     0.0360    0.0830    0.1785    0.4155    0.8385    0.9980    1.0000
%     0.0535    0.0935    0.1960    0.4605    0.9045    1.0000    1.0000
%     0.0585    0.1175    0.2025    0.4970    0.9225    1.0000    1.0000
%     0.0755    0.1145    0.2115    0.5000    0.9390    1.0000    1.0000
%     0.0750    0.1165    0.2285    0.5185    0.9445    1.0000    1.0000
%     0.0890    0.1205    0.2375    0.5235    0.9625    1.0000    1.0000

% method = 'nagao', reps = 2000; 25.11.14
% prob =
% 
%     0.0555    0.1050    0.1780    0.3995    0.7700    0.9795    1.0000
%     0.0815    0.0905    0.1950    0.4260    0.8470    0.9990    1.0000
%     0.0570    0.1005    0.1995    0.4605    0.9005    1.0000    1.0000
%     0.0435    0.1065    0.2150    0.4700    0.9335    1.0000    1.0000
%     0.0475    0.1035    0.2085    0.4815    0.9370    1.0000    1.0000
%     0.0535    0.1110    0.2175    0.5090    0.9490    1.0000    1.0000
%     0.0560    0.0995    0.1855    0.4890    0.9525    1.0000    1.0000

% method = 'wang', reps = 2000; 25.11.14
% prob =
% 
%     0.0580    0.1025    0.2190    0.5230    0.8735    0.9945    1.0000
%     0.1220    0.1930    0.2915    0.5815    0.9265    1.0000    1.0000
%     0.1375    0.2855    0.4100    0.6500    0.9625    1.0000    1.0000
%     0.1570    0.3140    0.4950    0.7425    0.9810    1.0000    1.0000
%     0.1565    0.3085    0.5345    0.8165    0.9900    1.0000    1.0000
%     0.1775    0.3290    0.5225    0.8010    0.9970    1.0000    1.0000
%     0.1560    0.3235    0.5290    0.8190    0.9965    1.0000    1.0000
% 
