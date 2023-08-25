%%%%%%%%%%%%%%% Numerical Examples: Comparing No-Borrowing and Leverage Economy %%%%%%%%%%%%%%%

%% 2 period, Leverage Economy 
%{
note that the system is nonlinear
syms p0 h
ic = h + (1-h)*(h+0.2*(1-h)) == p0;
mc = (1+p0)/p0*(1-h) == 1;

[A,b] = equationsToMatrix([ic, mc], [p0, h]);

x=linsolve(A,b)
%}

% h0= 0.6861; p= 0.7489
x = fsolve(@LeverageEconomy, [0.5,0.9]);
[h0,p] = struct('x', num2cell(x)).x;


%% 3-period, No Borrowing
% h= 0.545; p0= 0.835; pD= 0.636
x = fsolve(@NoBorrow, [1,1,1]);
[h0,p0,pD] = struct('x', num2cell(x)).x; crash=(p0-pD)/p0*100;
strcat("Price Crash in No Borrowing:",num2str(crash,4),"%")


%% 3-period, With Borrowing
% h0 = 0.8699, hD = 0.6165, p0 = 0.9465, pD = 0.6932
x = fsolve(@WithBorrow, [1,1,1,0.5]);
[h0,hD,p0,pD] = struct('x', num2cell(x)).x;crash=(p0-pD)/p0*100;
strcat("Price Crash in No Borrowing:",num2str(crash,4),"%") %>>> 26.76%
% try if no crash in fundamental
x = fsolve(@WithBorrow2, [0.7,0.6,0.9,0.5]);
[h0,hD,p0,pD] = struct('x', num2cell(x)).x;crash=(p0-pD)/p0*100;
strcat("Price Crash as if no belief crash:",num2str(crash,4),"%") %>>> 15.52%
strcat("Price Down as if no belief crash:",num2str(pD,4)) %>>> 15.52%


function F = WithBorrow2(x)
    % As if subjective belief of crash doesn't change (cannot fully explain
    % the price crash of 26.76%)

    R = 0.2;
    [h0,hD,p0,pD] = struct('x', num2cell(x)).x;
    qh0 = sqrt(1-h0); % qhD = hD;
    F(1) = (1-h0)*(1+p0) - (p0-pD);  % mc at 0 
    % opt for h0
    %                                                                                                                                           
    F(2) = qh0*(pD-R)*(p0-pD) + (1-qh0)*qh0*(1-R)*(p0-pD)-qh0*(1-pD)*(pD-R);   
    F(3) = (h0-hD)*(1+p0)-(pD-R);    % mc at D
    F(4) = hD + (1-hD)*R - pD;       % opt at D
end

function F = WithBorrow(x)
    R = 0.2;
    [h0,hD,p0,pD] = struct('x', num2cell(x)).x;
    F(1) = (1-h0)*(1+p0) - (p0-pD);  % MC at 0 
    % IC for h0
    F(2) = h0*(pD-R)*(p0-pD) + (1-h0)*h0*(1-R)*(p0-pD)-h0*(1-pD)*(pD-R);   
    F(3) = (h0-hD)*(1+p0)-(pD-R);    % mc at D
    F(4) = hD + (1-hD)*R - pD;       % opt at D
end


function F = NoBorrow(x)
    R=0.2;
    [h,p0,pD] = struct('x', num2cell(x)).x;

    F(1) = (1-h)*(1+p0) - p0;  % MC
    F(2) = h + (1-h)*h +(1-h)^2*R - p0;
    F(3) = h + (1-h)*R - pD;   % IC
end

function f = LeverageEconomy(x)
    R = 0.2;
    [h0,p] =  struct('x', num2cell(x)).x;

    % opt
    f(1) = h0*(1-R)-(p-R);

    % mc 
    f(2) = (1-h0)*(1+p) - (p-R);
end







