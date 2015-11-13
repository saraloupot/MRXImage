function header()
%HEADER
%file of constants that can be read in for use in various
%subroutines/macros

%% Exclude detectors 
exclude= [6];
%% Multiple stage positions
stgpos = [0 0 0; -1 0 0];
%% Dimensions
ndimx=25;
ndimy = 25;
ndimz=25;
xFOVLimits = [-3.5,3.5];
yFOVLimits = [-3.5,3.5];
zFOVLimits = [0,4.5];
xWorldLimits = [-5,5];
yWorldLimits = [-5,5];
zWorldLimits = [0,6];
xPhantomCoords = -3.6:.9:3.6;
yPhantomCoords = 3.6:-.9:-3.6;

xunits = linspace(xFOVLimits(1),xFOVLimits(2),ndimx);
yunits = linspace(yFOVLimits(1),yFOVLimits(2),ndimy);
zunits = linspace(zFOVLimits(1),zFOVLimits(2),ndimz);
xworldunits = linspace(xWorldLimits(1),xWorldLimits(2),ndimx);
yworldunits = linspace(yWorldLimits(1),yWorldLimits(2),ndimy);
zworldunits = linspace(zWorldLimits(1),zWorldLimits(2),ndimz);
[X, Y, Z] = ndgrid(xunits,yunits,zunits);
Xmgrid = permute(X,[2,1,3]);
Ymgrid = permute(Y,[2,1,3]);
Zmgrid = permute(Z,[2,1,3]);
points = [X(:),Y(:),Z(:)];
[Xworld, Yworld, Zworld] = ndgrid(xworldunits, yworldunits, zworldunits);
rdet = [-.5, -.4,5.5; -2.4,.6,5.5; -.6,1.6,5.5; 1.4, .9, 5.5; 1.2, -1.5, 5.5;-0.985,-1.706,5.5;  -2.4,-1.5,5.5];
%% exclude data
detectors = 1:7;
use = find(~ismember(detectors,exclude));
rdet = rdet(use,:);
%% add stage positions
rdet_temp = []; use_temp=[];
for ii = 1:size(stgpos,1)
    mvstage = repmat(stgpos(ii,:),size(rdet,1),1);
    rdet_temp = [rdet_temp; rdet-mvstage];
    use_temp = [use_temp, use+repmat((ii-1)*7,1,size(rdet,1))];
end
rdet=rdet_temp;
use =use_temp;
clear rdet_temp mvstage;
rdetfake = [rdet; -7, 7, rdet(1,3);7, 7, rdet(1,3); -7, -7, rdet(1,3); 7,-7,rdet(1,3) ];
% BFake = [B; 0; 0; 0; 0];
%% Calculate the unbiasing matrix

ndet=size(rdet,1);
% Calculate the biasing matrix
for detect = 1:ndet % Loop through the detectors
    for ii = 1:length(X(:))
        r_s = points(ii,:);
        rvec = (rdet(detect,:)-r_s)./100; % to convert to m
         W(detect,ii) = 1E-7*((3*rvec(1,3).^2)/norm(rvec,2)^5 - 1/norm(rvec,2)^3);
    end            
end % end detector loop
G = sqrt(sum(W.^2,1));
G1= diag(G.^-1);
A = W*G1;
clear G1;

%% save to a single variable
save('headerVariables');