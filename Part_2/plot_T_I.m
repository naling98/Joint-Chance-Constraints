clear all;
clc;
close all;

A = load('profile_1','T', 'I', 'mat_ne', 'all_payoffs');

B = load('profile_1_beta', 'T', 'I', 'mat_ne', 'all_payoffs');

A = struct2cell(A);
B = struct2cell(B);

n = size(A{2})(1,1);

Mean_T{1} = sum(A{2})/n;
Variance_T{1} = sum((A{2}-Mean_T{1}).^2)/n;
SD_T{1} = sqrt(Variance_T{1});

Mean_T{2} = sum(B{2})/n;
Variance_T{2} = sum((B{2}-Mean_T{2}).^2)/n;
SD_T{2} = sqrt(Variance_T{2});

disp(['Mean Time with Chance : ' num2str(Mean_T{1})]);
disp(['SD Time with Chance : ' num2str(SD_T{1})]);
disp(' ');
disp(['Mean Time without Chance : ' num2str(Mean_T{2})]);
disp(['SD Time without Chance : ' num2str(SD_T{2})]);
disp(' ');
disp('--------------------------------------');
disp(' ');

Mean_I{1} = sum(A{3})/n;
Variance_I{1} = sum((A{3}-Mean_I{1}).^2)/n;
SD_I{1} = sqrt(Variance_I{1});

Mean_I{2} = sum(B{3})/n;
Variance_I{2} = sum((B{3}-Mean_I{2}).^2)/n;
SD_I{2} = sqrt(Variance_I{2});

disp(['Mean Iterations with Chance: ' num2str(Mean_I{1})]);
disp(['SD Iterations with Chance: ' num2str(SD_I{1})]);
disp(' ');
disp(['Mean Iterations without Chance : ' num2str(Mean_I{2})]);
disp(['SD Iterations without Chance : ' num2str(SD_I{2})]);
disp(' ');
disp('--------------------------------------');
disp(' ');

s = size(A{4}{1})(1,2);

Er1 = zeros(s,1);
Er2 = zeros(s,1);

for i = 1:s
    C1 = A{4}{1}{i} - B{4}{1}{i};
    C2 = A{4}{2}{i} - B{4}{2}{i};

    Er1(i,1) = sqrt(mean(mean(C1.^2)));
    Er2(i,1) = sqrt(mean(mean(C2.^2)));
   
endfor

disp(['Player 1 av. NE difference between chance and non-chance : ' num2str(mean(Er1))]);
disp(['Player 2 av. NE difference between chance and non-chance : ' num2str(mean(Er2))]);
disp('');

z = 1:n;
plot(z,A{2})
hold
plot(z,B{2})
title('CHANGING INITIAL POINT in order');
xlabel('Instance');
ylabel('Time in sec');
legend("With Chance", "Without Chance");

figure;

plot(z,A{3})
hold
plot(z,B{3})
title('CHANGING INITIAL POINT in order');
xlabel('Instance');
ylabel('Iterations');
legend("With Chance", "Without Chance");

figure;

T = zeros(s,2);

for i=1:s
    T(i,1) = A{1}{i}{1}/1000;
    T(i,2) = A{1}{i}{2}/1000;
endfor

sd_per = (sqrt(sum((T - mean(T)).^2)/s)./mean(T));

disp(['Normalized SD in Payoffs : ' num2str(sd_per)]);
plot(z,T(:,1))
hold
plot(z,T(:,2))
title('Payoffs at Nash Eq. with CHANGING INITIAL POINT');
xlabel('Instance');
ylabel('Payoffs (in Thousands)');
legend("P1", "P2");
