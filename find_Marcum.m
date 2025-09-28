M = 3;      % M
a = 0;      % a
y_target = 10^(-4);  % Marcum Q value

fun = @(b) marcumq(0,max(0, b),M) - y_target;


b_initial = 7; % initial guess
b_solution = fzero(fun, b_initial);

disp(['solved b ', num2str(b_solution)]);

%%%  -4     5.2779 
%%%  -4.5   5.5224;
%%%  -5     5.7539;
%%% Pfa = 10^-5.5, b = 5.9744;
