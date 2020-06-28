clear all;
clc;
close all;

for t = 1:25
    A = load(['profile_1_instance' num2str(t)], 'payoffs_history');

    A = struct2cell(A);

    n = size(A{1}{1})(1,2);

    y = 1:min(n,20);
    subplot(5,5,t)
    plot(y,A{1}{1}(1:min(n,20))/1000);
    hold
    plot(y,A{1}{2}(1:min(n,20))/1000);
    title(['Instance ' num2str(t) ' with chance']);
    xlabel('Iteration');
    ylabel('Payoffs (in Thousands)');
 
endfor
figure;
for t = 1:25
    A = load(['profile_1_instance_beta' num2str(t)], 'payoffs_history');

    A = struct2cell(A);

    n = size(A{1}{1})(1,2);

    y = 1:n;
    subplot(5,5,t)
    plot(y,A{1}{1}/1000);
    hold
    plot(y,A{1}{2}/1000);
    title(['Instance ' num2str(t) ' without chance']);
    xlabel('Iteration');
    ylabel('Payoffs (in Thousands)');
 
endfor