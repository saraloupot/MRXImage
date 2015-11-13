function plotSourceScatter(s,B,varargin)
if(nargin==2)
    header;
     mkr='ro'; mkf=true;
elseif(nargin==3)
        X=varargin{1}.X;
        Y=varargin{1}.Y;
        Z=varargin{1}.Z;
        rdet=varargin{1}.rdet;
        xWorldLimits=varargin{1}.xWorldLimits;
        yWorldLimits=varargin{1}.yWorldLimits;
        zWorldLimits=varargin{1}.zWorldLimits;
        points=varargin{1}.points;
        zunits=varargin{1}.zunits;
        Xworld=varargin{1}.Xworld;
        Yworld=varargin{1}.Yworld;
        ndimx=varargin{1}.ndimx;
        ndimy=varargin{1}.ndimy;
        mkr='ro'; mkf=true;
else
    if(isempty(varargin{1}))
        header; 
    else
        X=varargin{1}.X;
        Y=varargin{1}.Y;
        Z=varargin{1}.Z;
        rdet=varargin{1}.rdet;
        xWorldLimits=varargin{1}.xWorldLimits;
        yWorldLimits=varargin{1}.yWorldLimits;
        zWorldLimits=varargin{1}.zWorldLimits;
        points=varargin{1}.points;
        zunits=varargin{1}.zunits;
        Xworld=varargin{1}.Xworld;
        Yworld=varargin{1}.Yworld;
        ndimx=varargin{1}.ndimx;
        ndimy=varargin{1}.ndimy;
    end
    if(isempty(varargin{2}))
        mkr='ro'; mkf=true;
    else
    mkr = varargin{2}.mkr; mkf=varargin{2}.mkf;
    end
end
Xcoords = s.source(:,1);
Ycoords = s.source(:,2);
Zcoords = s.source(:,3);
strength= s.source(:,4)./max(s.source(:,4));
rotate3d on; 
scatter3(rdet(:,1),rdet(:,2),rdet(:,3),200,B,'filled');rotate3d on; hold on;
if(mkf)
    scatter3(Xcoords,Ycoords,Zcoords,200.*strength,mkr,'filled');
else
    scatter3(Xcoords,Ycoords,Zcoords,200.*strength,mkr);
end
xlim([min(xWorldLimits) max(xWorldLimits)]); xlabel('x cm','fontsize',14);
zlim([0 max(zWorldLimits)]); ylabel('y cm','fontsize',14);
ylim([min(yWorldLimits) max(yWorldLimits)]); zlabel('z cm','fontsize',14);
view([0 0]); set(gca,'fontsize',18);
shading interp;
hold off;
