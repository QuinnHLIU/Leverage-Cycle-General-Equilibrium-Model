%%%%%%%%%%%%%%% Numerical Examples For Financial Innovation %%%%%%%%%%%%%%%



%% Leverage 
x = fsolve(@Leverage,[0.5,0.5]);
[h1,p] = struct('x', num2cell(x)).x
%{
h1 =
    0.6340
p =
    0.8928
%}

%% Tranching  
x = fsolve(@Tranching,[0.5,0.2,0.9,0.1]);
[h0,hT,p,piT] = struct('x', num2cell(x)).x
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

%% Financial Integration
R=0.2; x = fsolve(@FI,[0.61,0.04,1,0.9,0.1]);
[h0,hT,p,pstar,piT] = struct('x', num2cell(x)).x
%{
h0 =
    0.6096
hT =
    0.0465
p =
    1.0294
pstar =
    0.8780
piT =
    0.1818
%}
FbuyY = (1-h0)*(1+pstar)/(p-piT+pstar-R)*p; strcat("Foreign buy Y:",num2str(FbuyY,4))
FsellT= (1-h0)*(1+pstar)/(p-piT+pstar-R)*piT; strcat("Foreign sell Tranche:",num2str(FsellT,4))
HbuyYstar= (1-h0)*(1+p)/(p-piT+pstar-R)*pstar; strcat("Home buy Y*:",num2str(HbuyYstar,4))

% how many tranche Home sold to Foreign = how many Y* H bought from F
% how many debt Home bought from Foreign= how many Y  F bought from H
CAdeficit= FbuyY + (1-h0)*(1+p)/(p-piT+pstar-R)*piT...  % Home sell Y and Tranche (not in F in autarky)
    - HbuyYstar - (1-h0)*(1+pstar)/(p-piT+pstar-R)*R;   % Home buy Y* and Debt (not in H in autarky)
strcat("Home buy Y*:",num2str(CAdeficit,4))


%% Leverage Economy with 3 States
M=0.93; D=.81; 
x = fsolve(@Leverage3,[0.99,0.92,0.90,0.95,0.9]);
[hM,hD,hpi,p,piM] = struct('x', num2cell(x)).x;
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
%{
hM =
    0.9742
hpi =
    0.9231
p =
    0.9608 < hM
piM =
    0.9103                 
%}


function F = FI(x)
    
    R=0.2; 
    
    [h0,hT,p,pstar,piT] = struct('x', num2cell(x)).x;

    % opt of h0
    F(1) = gamma(h0) - p + piT;              
    F(2) = gamma(h0)*(1-R) - pstar + R;  

    % opt of hT
    F(3) = (1-gamma(hT))*R - piT;
    
    % market clearing
    F(4) = (1-h0)*(2+p+pstar) - p + piT - pstar + R;  % mc of Y,Y*
    F(5) = hT*(2+p+pstar) - piT;                      % mc of Tranche 

end

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

function g = gamma(h)
    g = 1-(1-h)^2;
end





