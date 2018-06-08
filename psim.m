function y = psim(mutRate, keyLen)
% calculate the probablity of the same key
% Lu Cheng
% 21.05.2018

y = (1-mutRate).^(2*keyLen);
