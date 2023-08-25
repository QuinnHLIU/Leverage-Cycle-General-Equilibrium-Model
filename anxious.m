%%%%%%%%% Dynamic Model with 3 periods and Multiple assets %%%%%%%%
%%%%%%%%% The model will showcase a spillover effect of bad news %%
%%% Further analysis on credit rating and primary issurance see Fostel and
%%% Geankoplos (2012). Note: two possibilities in the equlibrium regimes




%% Contagion with One EM Asset
% fsolve won't work fine if there are multiple equilibria
% but I checked with numerous ways, the solution fsolve finds is indeed the
% unique one
x_init = [0.9479,0.8506,0.8818,0.7055,0.7055];H=0.2; G=0.2; 
[p0H_init,p0G_init,pUG_init,pDH_init,pDG_init]=struct('x_init', num2cell(x_init)).x_init;
hD_init = (pDG_init-G)/(1-G) %>>>0.6319
hU_init = (pUG_init-G)/(1-G) %>>>0.8522
h0_init = hD_init + (pDH_init-H+pDG_init-G)/(1+p0H_init+p0G_init) %>>>0.9931
%h0_init = 1-(p0H-pDH+p0G-pDG)/(1+p0H+p0G) %>>>0.8266...
% then I found hU < h0 using values in Table 2.1; so it must be wrong
x = fsolve(@Contagion, [hU_init,h0_init,hD_init,p0H_init,p0G_init,pUG_init,pDH_init,pDG_init]);
[hU,h0,hD,p0H,p0G,pUG,pDH,pDG] = struct('x', num2cell(x)).x;
%0.8182    0.8266    0.5189    0.9036    0.7946    0.8545    0.6151    0.6151
crash=(p0H-pDH)/p0H*100;strcat("Price Crash of H:",num2str(crash,4),"%") %>>> 31.93%
crash=(p0G-pDG)/p0G*100;strcat("Price Crash of G:",num2str(crash,4),"%") %>>> 22.59%

% Also try Newton's method to solve systems of nonlinear algebraic
% equations, or I do grid search on initial guess. Whatever road I take, 
% I got exactly the same solution as fsolve. This just proves that the
% solution is indeed unique, and initial guess within reasonable region
% would just work fine.
%{
resultMat = zeros(50*50*50*50*50,8);
cnt=1;
for k1 = 1:20
    for k2 = 1:20
        for k3 = 1:20
            for k4 = 1:20
                for k5 = 1:20
p0H_init = k1*0.05; p0G_init = k2*0.05; pUG_init = k3*0.05;
pDH_init = k4*0.05; pDG_init = k5*0.05;
hD_init = (pDG_init-G)/(1-G);
hU_init = (pUG_init-G)/(1-G);
h0_init = hD_init + (pDH_init-H+pDG_init-G)/(1+p0H_init+p0G_init);
x = fsolve(@Contagion, [hU_init,h0_init,hD_init,p0H_init,p0G_init,pUG_init,pDH_init,pDG_init]);
[hU,h0,hD,p0H,p0G,pUG,pDH,pDG] = struct('x', num2cell(x)).x;
resultMat(cnt,:) = [hU,h0,hD,p0H,p0G,pUG,pDH,pDG];
cnt=cnt +1;
                end
            end
        end
    end
end

hUColumn = resultMat(:,1); h0Column = resultMat(:,2); hDColumn = resultMat(:,3);
p0HColumn = resultMat(:,4); p0GColumn = resultMat(:,5); pUGColumn = resultMat(:,6);
pDHColumn = resultMat(:,7); pDGColumn = resultMat(:,8);
rowsatisfy = unique(round(resultMat((p0HColumn>0)&(p0GColumn>0)&(pUGColumn>0)&(pDHColumn>0)&(pDGColumn>0), :),4),'rows');
% again, 0.8182	0.8266	0.5189	0.9036	0.7946	0.8545	0.6151	0.6151 is
% the only solution. This cannot be a consistent equilibrium because those
% who choose to hold X at 0 will have gross return = 1 instead of
% (1-pDH)/(p0H-pDH). The market clearing at U state would be wrong...and I
% tried to solve the equilibrium using the updated MC at U, the result is
% consistent with h0 > hU 
%}
[x, n] = SolveContagion_Newton(@Contagion, @JacobianContagion,...
        [hU_init,h0_init,hD_init,p0H_init,p0G_init,pUG_init,pDH_init,pDG_init], 0.0001);
[hU,h0,hD,p0H,p0G,pUG,pDH,pDG] = struct('x', num2cell(x)).x;
crash=(p0H-pDH)/p0H*100;strcat("Price Crash of H:",num2str(crash,4),"%")
crash=(p0G-pDG)/p0G*100;strcat("Price Crash of G:",num2str(crash,4),"%")


x_init = [0.9479,0.8506,0.8818,0.7055,0.7055];H=0.2; G=0.2; 
[p0H_init,p0G_init,pUG_init,pDH_init,pDG_init]=struct('x_init', num2cell(x_init)).x_init;
hD_init = (pDG_init-G)/(1-G); hU_init = (pUG_init-G)/(1-G);
h0_init = hD_init + (pDH_init-H+pDG_init-G)/(1+p0H_init+p0G_init) %>>>0.9931
%h0_init = 1-(p0H-pDH+p0G-pDG)/(1+p0H+p0G) %>>>0.8266...
%if hU turns out to be larger than h0, market clearing in U, AND OPT IN 0 must change !
x = fsolve(@Contagion2, [hU_init,h0_init,hD_init,p0H_init,p0G_init,pUG_init,pDH_init,pDG_init]);
[hU,h0,hD,p0H,p0G,pUG,pDH,pDG] = struct('x', num2cell(x)).x;
crash=(p0H-pDH)/p0H*100;strcat("Price Crash of H:",num2str(crash,4),"%") %>>> 31.83%
crash=(p0G-pDG)/p0G*100;strcat("Price Crash of G:",num2str(crash,4),"%") %>>> 22.72%

%{
hU =
    0.8218
h0 =
    0.8265
hD =
    0.5189
p0H =
    0.9023
p0G =
    0.7959
pUG =
    0.8574
pDH =
    0.6151
pDG =
    0.6151
%}



%% Anxious economy with Differential Contagion
% if hU > h0 >>> inconsistent
H=0.2; G=0.2; B=0.1;
p0H_i=0.9106;p0G_i=0.7786;p0B_i=0.7509; pUG_i=0.8226;pUB_i=0.8004;
pDH_i=0.6505;pDG_i=0.6505;pDB_i=0.6068;
h0_i = 1-(p0H_i-pDH_i+p0G_i-pDG_i+p0B_i-pDB_i)/(1+p0H_i+p0G_i+p0B_i);       %0.8453
hU_i = 1-(pUG_i-G+pUB_i-B)*(p0H_i-pDH_i)/(1+p0H_i+p0G_i+p0B_i)/(1-pDH_i);   %0.7138 < 0.8453 ???
hD_i = h0_i - (pDH_i-H+pDG_i-G+pDB_i-B)/(1+p0H_i+p0G_i+p0B_i);              %0.4360
x = fsolve(@diffContagion,[hU_i,h0_i,hD_i,p0H_i,p0G_i,p0B_i,pUG_i,pUB_i,pDH_i,pDG_i,pDB_i]);
[hU,h0,hD,p0H,p0G,p0B,pUG,pUB,pDH,pDG,pDB] = struct('x', num2cell(x)).x;
crash=(p0H-pDH)/p0H*100;strcat("Price Crash of H:",num2str(crash,4),"%") %>>> 35.13%
crash=(p0G-pDG)/p0G*100;strcat("Price Crash of G:",num2str(crash,4),"%") %>>> 21.12%
%{
hU =
    0.7278
h0 =
    0.8162
hD =
    0.4614
p0H =
    0.8773
p0G =
    0.7215
p0B =
    0.6521
pUG =
    0.7823
pUB =
    0.7550
pDH =
    0.5691
pDG =
    0.5691
pDB =
    0.5153
%}

% if hU < h0 >>> consistent
H=0.2; G=0.2; B=0.1;
p0H_i=0.9106;p0G_i=0.7786;p0B_i=0.7509; pUG_i=0.8226;pUB_i=0.8004;
pDH_i=0.6505;pDG_i=0.6505;pDB_i=0.6068;
h0_i = 1-(p0H_i-pDH_i+p0G_i-pDG_i+p0B_i-pDB_i)/(1+p0H_i+p0G_i+p0B_i);       %0.8453
hU_i = 1-(pUG_i-G+pUB_i-B)*(p0H_i-pDH_i)/(1+p0H_i+p0G_i+p0B_i)/(1-pDH_i);   %0.7138 < 0.8453 ???
hD_i = h0_i - (pDH_i-H+pDG_i-G+pDB_i-B)/(1+p0H_i+p0G_i+p0B_i);              %0.4360
x = fsolve(@diffContagion2,[hU_i,h0_i,hD_i,p0H_i,p0G_i,p0B_i,pUG_i,pUB_i,pDH_i,pDG_i,pDB_i]);
[hU,h0,hD,p0H,p0G,p0B,pUG,pUB,pDH,pDG,pDB] = struct('x', num2cell(x)).x
crash=(p0H-pDH)/p0H*100;strcat("Price Crash of H:",num2str(crash,4),"%") %>>> 34.07%
crash=(p0G-pDG)/p0G*100;strcat("Price Crash of G:",num2str(crash,4),"%") %>>> 22.49%
crash=(p0B-pDB)/p0B*100;strcat("Price Crash of B:",num2str(crash,4),"%") %>>> 22.78%
%{
hU =
    0.7633
h0 =
    0.8129
hD =
    0.4601
p0H =
    0.8616
p0G =
    0.7329
p0B =
    0.6657
pUG =
    0.8106
pUB =
    0.7869
pDH =
    0.5681
pDG =
    0.5681
pDB =
    0.5141
%}


%% Financial Innovation

%% Complete Markets
R=0.2;
pU_i=0.55; pD_i=0.45; pY_i=pU_i+R*pD_i; h1_i=1-pU_i*2/(1+pY_i);
x = fsolve(@CompleteMarket,[h1_i,pU_i,pD_i,pY_i]);
[h1,pU,pD,pY] = struct('x', num2cell(x)).x
%{
h1 =
    0.3292
pU =
    0.5500
pD =
    0.4500
pY =
    0.6400
%}


%% No Borrowing
x = fsolve(@NoBorrow2,[0.5,0.5]);
[h1,p] = struct('x', num2cell(x)).x;
%{
h1 =
    0.5451
p =
    0.8345
%}

%% Leverage 
x = fsolve(@Leverage,[0.5,0.5]);
[h1,p] = struct('x', num2cell(x)).x;
%{
h1 =
    0.6340
p =
    0.8928
%}

%% Tranching  
x = fsolve(@Tranching,[0.5,0.2,0.9,0.1]);
[h0,hT,p,piT] = struct('x', num2cell(x)).x;
%{
h0 =
    0.5852
hT =
    0.0841
p =
    0.9957
piT =
    0.1678
%}

%% CDS + Tranching  
x = fsolve(@CDS,[0.5,0.8,0.1,0.1]);
[h0,p,piT,piCDS] = struct('x', num2cell(x)).x
%{
exactly the complete market 
h0 =
    0.3292
p =
    0.6400
piT =
    0.0900
piCDS =
    0.3600
%}


%% Leverage Economy with 3 States
M=0.93; D=.81; 
x = fsolve(@Leverage3,[0.99,0.92,0.90,0.95,0.9]);
[hM,hD,hpi,p,piM] = struct('x', num2cell(x)).x;
iM=(M/piM-1)*100;strcat("Risk Premium:",num2str(iM,3),"%") 
RM=(p-piM)/p*100;strcat("Risky Margin:",num2str(RM,3),"%") 
SM=(p-D  )/p*100;strcat("Safe  Margin:",num2str(SM,3),"%") 
AM=((RM+SM)/2)  ;strcat("Avg Margin:"  ,num2str(AM,3),"%") 
%{
                                                                                                                                                                                         
hM =
    0.9951
hD =
    0.9567
hpi =
    0.8365
p =
    0.9093
piM =
    0.8734               
%}

%% Debt Collateralization
x = fsolve(@Pyramid,[0.99,0.92,0.90,0.95]);
[hM,hpi,p,piM] = struct('x', num2cell(x)).x;
iM=(M/piM-1)*100;strcat("Risk Premium:",num2str(iM,3),"%") 
RM=(p-piM)/p*100;strcat("Risky Margin:",num2str(RM,3),"%") 
SM=(p-D  )/p*100;strcat("Safe  Margin:",num2str(SM,3),"%") 
%{
hM =
    0.9742
hpi =
    0.9231
p =
    0.9608
piM =
    0.9103                 
%}




function F = Pyramid(x)

    M=0.93; D=.81; 
    
    [hM,hpi,p,piM] = struct('x', num2cell(x)).x;

    F(1) = gammaU(hM)*(1-M)*(piM-D) - (gammaU(hM)*(M-D)+gammaM(hM)*(M-D))*(p-piM); % opt of hM
    F(2) = gammaU(hpi)*(M-D)+gammaM(hpi)*(M-D) - (piM-D);                          % opt of hpi

    F(3) = (1-hM)*(1+p)-(p-piM);    % mc of Y
    F(4) = (hM-hpi)*(1+p)-(piM-D);  % mc of jM

end


function F = Leverage3(x)

    M=0.93; D=.81; 
    
    [hM,hD,hpi,p,piM] = struct('x', num2cell(x)).x;

    F(1) = gammaU(hM)*(1-M)*(p-D) - (gammaU(hM)*(1-D)+gammaM(hM)*(M-D))*(p-piM);           % opt of hM
    F(2) = gammaD(hD)*D+(1-gammaD(hD)*M)*(p-D) - (gammaU(hD)*(1-D)+gammaM(hD)*(M-D))*piM;  % opt of hD
    F(3) = gammaD(hpi)*D+(1-gammaD(hpi))*M - piM;        % opt of hpi

    F(4) = (1-hM)*(1+p)/(p-piM) + (hM-hD)*(1+p)/(p-D)-1; % mc of Y
    F(5) = (hD-hpi)*(1+p)*(p-piM) - (1-hM)*(1+p)*piM;    % mc of jM

end

function g = gammaD(h)
    g = 1-gammaU(h)-gammaM(h);
end
function g = gammaU(h)
    xi = 6.5;

    g = h^xi;
end
function g = gammaM(h)
    xi = 6.5;

    g = h^xi*(1-h^xi);
end

function F = CDS(x)
    R=0.2; 
    
    [h0,p,piT,piCDS] = struct('x', num2cell(x)).x;

    F(1) = R*piCDS - (1-R)*piT;            % indiff of D
    F(2) = 1-piCDS/(1-R) - p+piT;          % indiff of U

    F(3) = gamma(h0)*piT - (1-gamma(h0))*R*(p-piT);     %opt of h0 
    F(4) = (1-h0)*(1+p) - (p-piT+1-piCDS/(1-R)); %mc 

end

function F = Tranching(x)
    R=0.2; 
    
    [h0,hT,p,piT] = struct('x', num2cell(x)).x;

    F(1) = (1-h0)*(1+p) - p + piT;  % mc of Y 
    F(2) = hT*(1+p) - piT;          % mc of Tranche 

    F(3) = gamma(h0) - p + piT;     %opt of h0 
    F(4) = (1-gamma(hT))*R - piT;   %opt of hT

end


function F = Leverage(x)
    R=0.2; 
    
    [h1,p] = struct('x', num2cell(x)).x;

    F(1) = (1-h1)*(1+p) - p + R;
    F(2) = gamma(h1)*(1-R) - p + R;  %opt: lvg with minmax contract
end


function F = NoBorrow2(x)
    R=0.2;
    
    [h1,p] = struct('x', num2cell(x)).x;

    F(1) = (1-h1)*(1+p) - p;
    F(2) = gamma(h1)+ (1-gamma(h1))*R - p;
end



function F = CompleteMarket(x)
    R=0.2;
    
    [h1,pU,pD,pY] = struct('x', num2cell(x)).x;

    F(1) = (1-h1)*(1+pY) - 2*pU;
    F(2) = gamma(h1)*pD - (1-gamma(h1))*pU;
    F(3) = pY - pU - R*pD;
    F(4) = 1 - pU - pD;
end

function g = gamma(h)
    g = 1-(1-h)^2;
end


function F = diffContagion2(x)
    H=0.2; G=0.2; B=0.1;
    
    [hU,h0,hD,p0H,p0G,p0B,pUG,pUB,pDH,pDG,pDB] = struct('x', num2cell(x)).x;
    % equalized return
    F(1) = (1-pDH)*(p0G-pDG) - (pUG-pDG)*(p0H-pDH); % H and G at 0 
    F(2) = (pUG-pDB)*(p0B-pDB)-(pUB-pDB)*(p0G-pDG); % G and B at 0
    F(3) = (1-H)*(pDG-G) - (pDH-H)*(1-G);           % H and G at D 
    F(4) = (1-G)*(pDB-B)-(1-B)*(pDG-G);             % G and B at D 
    F(5) = (1-G)*(pUB-B)-(1-B)*(pUG-G);             % G and B at U 

    % market clearing
    F(6) = (1-h0)*(1+p0H+p0G+p0B)-(p0H-pDH+p0G-pDG+p0B-pDB);  % at 0 
    F(7) = (h0-hD)*(1+p0H+p0G+p0B)-(pDH-H+pDG-G+pDB-B);       % at D 
    F(8) = (1-hU)*(1+p0H+p0G+p0B)*(1-pDH)+(h0-hU)*(1+p0H+p0G+p0B)*(p0H-pDH)-(pUG-G+pUB-B)*(p0H-pDH);  % at U
    
    % indifference cond
    F(9) = hD*(1-H)-(pDH-H);  % at D 
    F(10)= hU*(1-G)-(pUG-G);  % at U
    F(11)= h0*(1-pDH)*(pDH-H)*(pUG-G)-h0^2*(1-G)*(pDH-H)*(p0H-pDH)-(1-h0)*h0*(1-H)*(p0H-pDH)*(pUG-G); % at 0
end


function F = diffContagion(x)
    H=0.2; G=0.2; B=0.1;
    
    [hU,h0,hD,p0H,p0G,p0B,pUG,pUB,pDH,pDG,pDB] = struct('x', num2cell(x)).x;
    
    F(1) = (1-pDH)*(p0G-pDG) - (pUG-pDG)*(p0H-pDH);
    F(2) = (1-H)*(pDG-G) - (pDH-H)*(1-G);
    F(3) = (1-h0)*(1+p0H+p0G+p0B)-(p0H-pDH+p0G-pDG+p0B-pDB);
    F(4) = (h0-hD)*(1+p0H+p0G+p0B)-(pDH-H+pDG-G+pDB-B);
    F(5) = (1-hU)*(1+p0H+p0G+p0B)*(1-pDH)-(pUG-G+pUB-B)*(p0H-pDH);
    F(6) = hD*(1-H)-(pDH-H);
    F(7) = hU*(1-G)-(pUG-G);
    F(8) = h0*(1-pDH)*(pDH-H)-h0*(pDH-H)*(p0H-pDH)-(1-h0)*h0*(1-H)*(p0H-pDH);
    F(9) = (pUG-pDB)*(p0B-pDB)-(pUB-pDB)*(p0G-pDG);
    F(10)= (1-G)*(pDB-B)-(1-B)*(pDG-G);
    F(11)= (1-G)*(pUB-B)-(1-B)*(pUG-G);

end

function J = JacobianContagion(x)
    % first derivative w.r.t. x(i) on the ith column
    H=0.2; G=0.2; 
    [hU,h0,hD,p0H,p0G,pUG,pDH,pDG] = struct('x', num2cell(x)).x;
    % use sympy in python to compute Jacobian (see calJacob.py)
    %{
        >>> from sympy import *
        >>> x0, x1 = symbols('x0 x1')
        >>> F0 = x0**2 - x1 + x0*cos(pi*x0)
        >>> F1 = x0*x1 + exp(-x1) - x0**(-1)
        >>> diff(F0, x0)
    %}
    J = [0,0,0,pDG - pUG,1 - pDH,-p0H + pDH,-p0G + pUG,p0H - 1;...
    0,0,0,0,0,0,G - 1,1 - H;...
    0,-p0G - p0H - 1,0,-h0,-h0,0,1,1;...
    0,p0G + p0H + 1,-p0G - p0H - 1,h0 - hD,h0 - hD,0,-1,-1;...
    -(1 - pDH)*(p0G + p0H + 1),0,0,G - pUG + (1 - hU)*(1 - pDH),(1 - hU)*(1 - pDH),-p0H + pDH,-G + pUG - (1 - hU)*(p0G + p0H + 1),0;...
    0,0,1 - G,0,0,0,0,-1;...
    1 - G,0,0,0,0,-1,0,0;...
    0,h0*(1 - H)*(p0H - pDH) - (1 - H)*(1 - h0)*(p0H - pDH) + (1 - pDH)*(-H + pDH) - (-H + pDH)*(p0H - pDH),...
    0,-h0*(1 - H)*(1 - h0) - h0*(-H + pDH),0,0,h0*(1 - H)*(1 - h0) + h0*(1 - pDH) - h0*(p0H - pDH),0];
    %J = sum(J,2);  %??
end


function F = Contagion2(x)
    % Rewrite the equations as if h0 > hU 

    H=0.2; G=0.2; 
    [hU,h0,hD,p0H,p0G, pUG,pDH,pDG] = struct('x', num2cell(x)).x;  %pUH...but there are many 1s 
    
    % equalize return
    F(1) = (1-pDH)/(p0H-pDH) - (pUG-pDG)/(p0G-pDG);  % at 0 
    F(2) = (1-H)/(pDH-H) - (1-G)/(pDG-G);            % at D
    
    % market clearing
    F(3) = (1 -h0)*(1+ p0H + p0G) -(p0H-pDH+p0G-pDG);         % at 0
    F(4) = (h0-hD)*(1+ p0H + p0G) -(pDH- H +pDG-  G);         % at D
    %F(5) = (1 -hU)*(1+ p0H + p0G)*(1-pDH)-(pUG-G)*(p0H-pDH); % at U
    F(5) = (1 -hU)*(1+ p0H + p0G)*(1-pDH)+(h0-hU)*(1+ p0H + p0G)*(p0H-pDH)-(pUG-G)*(p0H-pDH);

    % indifference condition
    F(6) = hD*(1-G) -(pDG-G);
    F(7) = hU*(1-G) -(pUG-G); 
    %F(8) = h0*(1-pDH)*(pDH-H)-h0*(pDH-H)*(p0H-pDH)-(1-h0)*h0*(1-H)*(p0H-pDH);
    F(8) = h0*(1-pDH)*(pDH-H)*(pUG-G)-h0^2*(1-G)*(pDH-H)*(p0H-pDH)-(1-h0)*h0*(1-H)*(p0H-pDH)*(pUG-G);
end


function F = Contagion(x)
    H=0.2; G=0.2; 
    [hU,h0,hD,p0H,p0G, pUG,pDH,pDG] = struct('x', num2cell(x)).x;  %pUH...but there are many 1s 
    
    % equalize return
    F(1) = (1-pDH)*(p0G-pDG) - (pUG-pDG)*(p0H-pDH);  % at 0 
    F(2) = (1-H)*(pDG-G) - (1-G)*(pDH-H);            % at D
    
    % market clearing
    F(3) = (1 -h0)*(1+ p0H + p0G) -(p0H-pDH+p0G-pDG);        % at 0
    F(4) = (h0-hD)*(1+ p0H + p0G) -(pDH- H +pDG-  G);        % at D
    F(5) = (1 -hU)*(1+ p0H + p0G)*(1-pDH)-(pUG-G)*(p0H-pDH); % at U

    % indifference condition
    F(6) = hD*(1-G) -(pDG-G);  % at D 
    F(7) = hU*(1-G) -(pUG-G);  % at U 
    F(8) = h0*(1-pDH)*(pDH-H)-h0*1*(pDH-H)*(p0H-pDH)-(1-h0)*h0*(1-H)*(p0H-pDH);  % at 0
end

function [x, iteration_counter] = SolveContagion_Newton(F, J, x, eps)
    % Solve nonlinear system F=0 by Newton's method.
    % J is the Jacobian of F. Both F and J must be functions of x.
    % At input, x holds the start value. The iteration continues
    % until ||F|| < eps.
    
    F_value = F(x);
    F_norm = norm(F_value);  % l2 norm of vector
    iteration_counter = 0;
    while abs(F_norm) > eps && iteration_counter < 100
        delta = J(x)\(-F_value');  %=inv(J)*(-F') >>> n x 1 
        x = (x' + delta)';         %>>> 1 x n 
        F_value = F(x);            %>>> 1 x n
        F_norm = norm(F_value);
        iteration_counter = iteration_counter + 1;
    end

    % Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps
        iteration_counter = -1;
    end

end


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







