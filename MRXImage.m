function [varargout]=MRXImage(varargin)
%% MRXImage is a function to reconstruct a source distribution image from NanoMRX data
%% initialization
% Optinal name,value pairs: [default]
%   1. Source: [] if empty, run makeBField.  Otherwise, provide a structure
%   with B
%       fields: sourcelocation, sourcestrength, B and truth;
%   2. Algorithms: ['Linprog'],'sparseApproxLinprog','FOCUSS', 'CVX_penalty',
%       'CVX_error', 'CVX_Focuss'
%   3. Plotting: ['off'],'line', 'contour', 'scatter', 'montage', 'stats'
%
% Parameters: Algorithm-specific parameters include:
%       L: [1000] lambda value for FOCUSS and CVX_penalty
%       P: [0.1] p norm for FOCUSS
%       tol: [0.01*norm(B)] residual tolerance for CVX_error and CVX_Focuss
%       TO DO: zthresh (cvx_Delta): [0.01*norm(B)] threshold for zero in CVX_Focuss
% Options:
%       Output: [false] print results to cmd window
%       Save: ['filename'] save files to 'filename'
%       DrawContour: [false] draw contours of A on scatter plots
%       DrawSolution: [false] plot the true solution
%       DrawPhantom: [false] draw the phantom on scatter plots

% Uses the header and makeBField that is in the folder containing MRXImage.m
dstr=datestr(now,'yy_mm_dd_HHMM');
%Parse inputs
input=inputParser;
% Defaults:
Default_source = [];
Default_tol = [];
Default_DrawPhantom = false;
Default_DrawSolution = false;
Default_DrawContour = false;
Default_Output = false;
Default_maxIter = 15;
Default_plotname=['MRXImageResults_',dstr];
Default_Save = [];
Default_SNR = [];
Default_truth = [];
Default_quiet = true;
Default_headerVariables= 'headerVariables';
% TO DO: checkSave= @(x) any([isempty(x), ischar(x)]);
Default_Algorithm = {'SeDuMi'};
validAlgorithms = {'SeDuMi','SDPT3'};
% TO DO: checkAlgorithm = @(x) any(validatestring(x,validAlgorithms));

Default_plotting = {'off'};
validPlotting = {'off','line','contour','scatter','montage','stats'};
% TO DO: checkPlotting = @(x) any(validatestring(x,validPlotting));

% Add options to p
% TO DO: make first three addptional instead of parameters
addParameter(input,'B',Default_source);
addParameter(input,'Algorithm',Default_Algorithm);
addParameter(input,'Plotting',Default_plotting);
addParameter(input,'tol',Default_tol,@isnumeric);
addParameter(input,'DrawPhantom',Default_DrawPhantom,@islogical);
addParameter(input,'DrawContour',Default_DrawContour,@islogical);
addParameter(input,'Output',Default_Output,@islogical);
addParameter(input,'Save',Default_Save);
addParameter(input,'DrawSolution',Default_DrawSolution,@islogical);
addParameter(input,'Plotname',Default_plotname);
addParameter(input,'SNR',Default_SNR);
addParameter(input,'truth',Default_truth);
addParameter(input,'maxIter',Default_maxIter);
addParameter(input,'quiet',Default_quiet);
addParameter(input,'headerVariables',Default_headerVariables);
% Parse input and assign values
parse(input,varargin{:});
headerFile=input.Results.headerVariables;
load(headerFile);
B = input.Results.B;
truth=input.Results.truth;
if(isempty(truth))
    sourcelocation = [0 0 0];
    sourcestrength = [0 0 0];
    existsTruth=false;
else
    sourcelocation = truth.sourcelocation;
    sourcestrength = truth.sourcestrength;
    existsTruth=true;
end
if(isempty(B)); makeBField; end
Algorithms = input.Results.Algorithm;
plotting = input.Results.Plotting;
MRXsnr = input.Results.SNR;
maxIter = input.Results.maxIter;
cvx_tol = input.Results.tol;
if(isempty(cvx_tol)); cvx_tol=0.01*norm(B); end
drawPhantomOpt = input.Results.DrawPhantom;
drawContourOpt = input.Results.DrawContour;
drawSolution = input.Results.DrawSolution;
plotname=input.Results.Plotname;
output = input.Results.Output;
quiet = input.Results.quiet;
% TO DO: savefile = input.Results.Save;

% Print the true solution
if(output)
    clc;
    if(exist(truth,'var'))
    fprintf('%s\n',truth.out);
    else
        fprintf('There is no truth :( \n');
    end
end

%% Linear solves
for ialg=1:length(Algorithms)
    switch Algorithms{ialg}
      
        case 'SeDuMi'
            cvx_solver sedumi
            disp('SeDuMi');            
        case 'SPDT3'
            cvx_solver sdpt3
            disp('SDPT3');
    end
            Results=[];
            n=size(G,2);
            U=ones(n,1);
            cvx_optval_old = 1E10;
            Results=repmat(struct('name',[],'residualAy',[],'x_log',[],'x',[],'snr',[],'source',[],'id',[],'geodistance',[],'optval',[],'B',B,'tol',cvx_tol,'out',[],'U',[]),maxIter,1);
            for k= 1:maxIter
                if(quiet)
                    cvx_begin quiet
                else
                    cvx_begin
                end
                variable x_log(n)
                minimize(sum(U.*abs(x_log)))
                subject to
                norm(A*x_log-B,2) <=cvx_tol;
                x_log>=0;
                cvx_end
                cvxResults={cvx_optval; cvx_optbnd; cvx_slvitr; cvx_slvtol; cvx_status; cvx_cputime};
                cvx_Delta=0.01*max(x_log);
                U = 1./(cvx_Delta + abs(x_log));
                Results(k).residualAy=norm(A*x_log-B,2);
                Results(k).x_log=x_log;
                x_log=x_log./G';
                Results(k).x=x_log;
                Results(k).name=Algorithms{ialg};
                id=find(x_log>0.01*max(x_log));
                Results(k).id = id;
                Results(k).snr = MRXsnr;
                Results(k).source=[points(id,1),points(id,2),points(id,3),x_log(id)];
                if(sum(sourcestrength(:,3))>0)
                    Results(k).geodistance=gaussianDistance(Results(k).source,sourcelocation,sourcestrength);
                else
                    Results(k).geodistance = [];
                end
                Results(k).optval=cvx_optval;
                Results(k).out{1}=sprintf('===================================');
                Results(k).out{2}=sprintf('             %s iteration %i',Results(k).name,k);
                Results(k).out{3}=sprintf('===================================');
                Results(k).out{4}=sprintf('tol:                      %.2e', cvx_tol);
                Results(k).out{5}=sprintf('Solver Status               %s',cvx_status);
                Results(k).out{6}=sprintf('Residual ||Ax-b||:        %.2e',Results(k).residualAy);
                Results(k).out{7}=sprintf('SNR:                     %i',MRXsnr);
                Results(k).out{8}=sprintf('Geometric distance:    %10.2e\n',Results(k).geodistance');
                Results(k).out{9}=sprintf('===================================');
                Results(k).out{10}=sprintf('      X        Y       Z        mu ');
                Results(k).out{11}=sprintf('%8.2f%8.2f%8.2f%12.2e\n',Results(k).source');
                Results(k).U=U;
                if(abs(cvx_optval_old-cvx_optval)<0.001 || isnan(cvx_optval))
                    break;
                end
                cvx_optval_old = cvx_optval;
            end
            
            switch Algorithms{ialg}
                case 'SeDuMi'
                    SeDuMi_Results=Results(1:k);
                    assignin('caller','SeDuMi_Results',Results(1:k));
                case 'SDPT3'
                    SDPT3_Results=Results(1:k);
                    assignin('caller','SDPT3_Results',Results(1:k));
            end
            clear Results;
end
%% plotting

if(~strcmp(plotting,'off'))
    numCols = length(Algorithms);
    numRows = length(plotting);
    figure('Position',[500 500,1500,1000]);
    for iplotCol = 1:numCols
        alg = Algorithms{iplotCol};
        switch alg
            case 'SeDuMi'
                s = SeDuMi_Results(end);
            case 'SDPT3'
                s = SDPT3_Results(end);
        end
        subplot(numRows,numCols,iplotCol);
        for iplotRow = 1:numRows
            plottype = plotting{iplotRow};
            iplot = (iplotRow-1)*numCols+iplotCol;
            subplot(numRows,numCols,iplot);
            switch plottype
                case 'scatter'
                    hold on;
                    if(drawPhantomOpt)
                        drawPhantom;
                    end
                    if(drawContourOpt)
                        A1=reshape(A(1,:),ndimx,ndimy,ndimz);
                        A1mgrid = permute(A1,[2,1,3]);
                        hcont = contourslice(Xmgrid,Ymgrid,Zmgrid,A1mgrid.*B(1),[],sourcelocation(1,2),[]);
                        colormap(jet); caxis(caxis);
                    end
                    if(drawSolution)
                        scatter3(sourcelocation(:,1),sourcelocation(:,2), sourcelocation(:,3),250*sourcestrength/max(sourcestrength),'sk','LineWidth',1.5);
                    end
                    hold on;
                    marker.mkr='ro'; marker.mkf=true;
                    plotSourceScatter(s,B,headerVariables,marker);
                    hs=findobj(gca,'type','scatter');
                    hleg=legend([hs(3),hs(1)],'True location','Reconstructed location','location','northeast');
                    set(hleg,'fontsize',18);
                    
                case 'line'
                    plot(s.x./max(s.x),'-o'); hold on;
                    if(existsTruth) 
                        plot(truth.x./max(truth.x)); 
                        hleg=legend(s.name,'truth');
                        set(hleg,'fontsize',18);
                    end
                case 'stats'
                    if(existsTruth)
                        txt= text(0,0,[truth.out(:); s.out(:)],'fontsize',18);
                    else
                        txt= text(0,0,s.out(:),'fontsize',18);
                    end
                    set(gca,'visible','off');
                    xlim([0,1]); ylim([-1 1]);
            end
        end
        subplot(numRows,numCols,iplotCol);
        title(s.name,'fontsize',24);
    end
end

%% output
if(output)
    fprintf('=========================\n       End results: \n=========================\n');
    for ialg = 1:length(Algorithms)
        alg = Algorithms{ialg};
        switch alg
            case 'SeDuMi'
                s = SeDuMi_Results(end);
            case 'SDPT3'
                s = SDPT3_Results(end);
        end
        fprintf(s.out(:));
            fprintf('=====================================\n\n');
    end
    fprintf('=========================\n       All results: \n=========================\n');
    for ialg = 1:length(Algorithms)
        alg = Algorithms{ialg};
        switch alg
            case 'SeDuMi'
                s = SeDuMi_Results;
            case 'SDPT3'
                s = SDPT3_Results;
        end
        for iresult = 1:length(s)
            fprintf(s(iresults).out(:));
            fprintf('=====================================\n\n');
        end
    end
end
