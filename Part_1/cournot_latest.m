clear all;
clc;
close all;

warning('off','all');

pkg load statistics
pkg load optim


disp("Two firms in a Cournot Competetion with and without Chance Constraints");
disp(' ');

%-------------------------------------------------------------------------------------%
% Game Setup


global N N1 N2 N3 M x y z Ni A_u A_v A_s b x_min x_max beta delta c alpha mat_ne chance

for chance = 0:1
	
	S = 20; % Total number of instances
	T = zeros(S,1); %To store time for different data sizes
	I = zeros(S,1);
    all_payoffs = {};
	count = 1;
	mat_ne = {};
    
    N1 = 0; %N1: Solo plants for P1 
    N2 = 0; %N2: Solo plants for P2
    N3 = 4; %N3: Colocated Plants
    N = N1+N2+N3;  %N is the number of generation nodes
    Ni = {N1+N3, N2+N3};
    
    disp('................................................................');
    disp(['Number of Generation Nodes : ' num2str(N) ' [ Solo P1: ' num2str(N1) ', Solo P2: ' num2str(N2) ', Shared: ' num2str(N3) ']']);

    %N1+N2+N3 = N
    %1        - N1       : P1
    %N1+1     - N1+N3    : P1+P2
    %N1+N3+1  - N1+N2+N3 : P2

    M = 3; %M is the number of distribution nodes

    disp(['Numeber of Distribution Nodes : ' num2str(M)]);
    disp('................................................................');   
    disp(' ');
	
	for s = 1:S
		disp(['Chance = ' num2str(chance) ' | Instance : ' num2str(s) '  . . . ' ]);
		disp(' ');
		rand("seed", s);


		%-------------------------------------------------------------------------------------%
		% Production Setup
        
        	beta = 300 * ones(N,M);  %B = 100
		delta = 1 * ones(N,M); %δ = 1

		c = {0,0}; % Marginal Cost
		c{1} = 20 * ones(N1+N3, M);
		c{2} = 15 * ones(N2+N3, M);

		x_min = {0,0};
		x_min{1} = zeros(N1+N3, M);
		x_min{2} = zeros(N2+N3, M);

		bound1 = (beta(N1+1:N1+N3,:)-c{1}(N1+1:N1+N3,:))./delta(N1+1:N1+N3,:);
        	bound2 = (beta(N1+1:N1+N3,:)-c{2}(1:N3,:))./delta(N1+1:N1+N3,:);
        
		x_min{1}(1:N1,:) = (beta(1:N1,:)-c{1}(1:N1,:))/4;
		x_min{1}(N1+1:N1+N3,:) = (4*bound1 - bound2)/15;

		x_min{2}(N3+1:N2+N3,:) = (beta(N1+N3+1:N,:)-c{2}(N3+1:N2+N3,:))/4;
		x_min{2}(1:N3,:) = (4*bound2 - bound1)/15;

		x_max = {0,0};
		x_max{1} = 200 + 10*rand(N1+N3, M);
		x_max{2} = 200 + 10*rand(N2+N3, M);
		

##        disp('Production Range for P1 : ');
##        disp([x_min{1} x_max{1}]);
##        disp(' ');
##        disp('Production Range for P2 : ');
##        disp([x_min{2} x_max{2}]);
##        disp(' ');


		%-------------------------------------------------------------------------------------%
		% Distribution Setup

		x = {0,0};
		x{1} = 12*ones(Ni{1}, M);
		x{2} = 12*ones(Ni{2}, M);

		
		y = {log(x{1}), log(x{2})};
		z = {ones(Ni{1},1)/Ni{1}, ones(Ni{2},1)/Ni{2}};

		alpha = {0.9, 0.9};

		% There are Ni{i} chance constraints for each player
        if chance == 1
            A_u = {};
            A_v = {};  %A_v is the convariance matrix, positive semidefinite
            A_s = {};  %A_s is the square root of the above, but useful for calc., and used in generation anyway
            b = {};
            
            for i = 1:2
            A_u{i} = 0.2 + rand(Ni{i}, M)/10;
                b{i} = 150 + 40 * rand(Ni{i}, 1);
                for j = 1:Ni{i}
                    A_s{i}{j} = rand(M, M)/5;
                    A_v{i}{j} = A_s{i}{j}*A_s{i}{j}';
                endfor
            endfor 
        endif   

		%-------------------------------------------------------------------------------------%
		% Payoff et Utility

		function vals = utility(x, printPrice = false)
			vals = {0,0};

			global N N1 N2 N3 M beta delta c;

			supply_1 = [x{1}; zeros(N2,M)]; 
			supply_2 = [zeros(N1,M); x{2}];
			supply = supply_1 + supply_2;
			price = beta - delta.*supply;
            
			vals{1} =  sum(sum((price(1:N3+N1,:) - c{1}).*x{1}));
			vals{2} =  sum(sum((price(N1+1: N,:) - c{2}).*x{2}));            
		endfunction

        
		%-------------------------------------------------------------------------------------%
		% Vectorization Helpers

		w = [[reshape(y{1}, Ni{1}*M, 1); z{1}]; [reshape(y{2}, Ni{2}*M, 1); z{2}]] ;

		function z = extract_z(w)
			global N Ni M;
			z = {};
			z{1} = w((Ni{1}*M) +1 :(Ni{1}*M) + Ni{1});
			z{2} = w((Ni{1} + Ni{2})*M + Ni{1}+1 : end);
		endfunction

		function y = extract_y(w)
			global N Ni M;
			y = {};
			y{1} = reshape(w(1:Ni{1}*M), Ni{1}, M);
			y{2} = reshape(w(Ni{1}*M + Ni{1} + 1  :  (Ni{1} + Ni{2})*M + Ni{1}), Ni{2}, M);
		endfunction

		function x = extract_x(w)
			y = extract_y(w);
			x = {exp(y{1}), exp(y{2})};
		endfunction

		function w = create_w(wi, i)

			global y z Ni M;
			w = wi;
			j = 2/i;

			if i == 1
		  		w = [wi; reshape(y{j}, Ni{j}*M, 1) ; z{j}];
			elseif i == 2
		  		w = [reshape(y{j}, Ni{j}*M, 1); z{j}; wi];
			end
		endfunction
		%-------------------------------------------------------------------------------------%
		% Nonlinear Optimization
		function val = opti_payoff(w)
			val = utility(extract_x(w));
		endfunction


		function val = opti_g(w, i)
			z = extract_z(w){i};
			val = [1 - sum(z)];
		endfunction

		function val = opti_h(w, i)
			global M Ni x_min x_max b A_u A_s alpha chance;
			val = [];

			y = extract_y(w){i};
			z = extract_z(w){i};
			x = extract_x(w){i};

			val = [val; z];
			val = [val; reshape(x - x_min{i}, Ni{i}*M,1)];
            val = [val; reshape(x_max{i} - x, Ni{i}*M,1)];
			
            if chance == 1
                for k = 1:Ni{i}
                    t = b{i}(k,1);
                    t -= (A_u{i}(k,:))*x(k, :)';                    
                    p = alpha{i}**z(k,1);
                    
                    p = log(stdnormal_inv(p)); #Normal Dist Chance 
                    
                    p = p*ones(M, 1) + y(k,:)';
                    p = norm(A_s{i}{k} * exp(p));
                    t -= p;
                    val = [val; t];
                endfor
            endif  
		endfunction



		w = create_w([reshape(y{1}, Ni{1}*M, 1); z{1}], 1);
		starting_payoffs = opti_payoff(w);
		%-------------------------------------------------------------------------------------%

		maxIter = 5000;
		cool = true;
		tol = 1e-4;

		payoffs_history = {[starting_payoffs{1}], [starting_payoffs{2}]};

		tic; %starting time

		total_iter = 0;

		for iteration = 1:maxIter
			if (!cool)
				break
			endif

			for i = 1:2

                temp = x;
                vec = [reshape(y{i}, Ni{i}*M, 1); z{i}];
                [abc, obj, info, iter, nf, lambda] = sqp(vec,@(vec) -opti_payoff(create_w(vec, i)){i}, @(vec) opti_g(create_w(vec, i), i), @(vec) opti_h(create_w(vec, i), i), -realmax, +realmax);

                w = create_w(abc,i);
                x = extract_x(w);
                y = extract_y(w);
                z = extract_z(w);


                payoffs = opti_payoff(w);
                payoffs_history{1} = [payoffs_history{1}, payoffs{1}];
                payoffs_history{2} = [payoffs_history{2}, payoffs{2}];
			  


                disp(['Iteration Number : ' num2str(2*iteration+i-2)]);
                disp(['iter : ' num2str(iter)]);
                disp(mean(mean(abs(x{1}-temp{1}))));
                disp(mean(mean(abs(x{2}-temp{2}))));
                disp('----------------------');

                          

				if mean(mean(abs(x{1}-temp{1})))<tol && mean(mean(abs(x{2}-temp{2})))<tol
                                   
					opti = {[opti_g(w,1);opti_h(w,1)], [opti_g(w,2);opti_h(w,2)]};

					 if(~(-0.01<opti{1} && -0.01<opti{2}))
						total_iter = 2*iteration+i-2;
						disp('NO FEASIBLE SOLUTION');
						disp(' ');
						disp(['Total Number of Iterations : ' num2str(total_iter)]);  
                    	disp(' ');
    
					 else
						total_iter = 2*iteration+i-2;
						disp(['Total Number of Iterations : ' num2str(total_iter)]);
                        disp(' ');

					endif
					mat_ne{1}{s} = x{1};
					mat_ne{2}{s} = x{2}; 
					cool = false;
					break;
				endif
			endfor
			total_iter = 2*iteration+i-2;
		endfor

        all_payoffs{count} = payoffs;
		T(count) = toc;
		I(count) = total_iter;
		count++;
        disp('----------------------------------------------------------------------');
        disp('');

        if chance == 1
            %When Chance contraints are used
            save(['profile_1_instance' num2str(s)], 'N1','N2','N3','M','beta','delta','c','A_u','A_v','A_s','b','alpha','x_min','x_max','x','w','payoffs_history','total_iter', 'opti');
        else
            save(['profile_1_instance_beta' num2str(s)], 'N1','N2','N3','M','beta','delta','c','A_u','A_v','A_s','b','alpha','x_min','x_max','x','w','payoffs_history','total_iter', 'opti');
        endif
	endfor
	
	if chance == 1    
        %When Chance contraints are used
		save('profile_1','N1','N2','N3','M', 'all_payoffs', 'T','I', 'mat_ne');
	else
		save('profile_1_beta','N1','N2','N3','M', 'all_payoffs', 'T','I', 'mat_ne');
	endif
	
endfor
