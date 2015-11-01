function [original,npts]=analyze_everydata()
%
% analyze_data 
% read stat_datax, which contains all errors data and their matching
% parameters. And sort the errors according to their parameters which
% themselves have been binned into descrete numbers.
%
% OUTPUT: the averages and standard deviations for different bins of parameters.
% In the first step the parameters are the square of cosine and the cost
% function. The are both divided into 10 bins.
% Programmed by: Huagang Yan, yanhg@ccmu.edu.cn
% Aug 10, 2010
%--------------------------------------------------------------------------

% the upperlimit of npoints is 400,000,

original=zeros(1200000,10); % since each folder may contain different numbers of data points
                           % here we set an upper limit.
sorted1=zeros(10,10,30000);
sorted2=zeros(10,10,30000);
sorted3=zeros(10,10,30000);
sorted4=zeros(10,10,30000);
sorted5=zeros(10,10,30000);
sorted6=zeros(10,10,30000);
sorted7=zeros(10,10,30000);
sorted8=zeros(10,10,30000);
aves1=zeros(10,10);
aves2=zeros(10,10);
aves3=zeros(10,10);
aves4=zeros(10,10);
aves5=zeros(10,10);
aves6=zeros(10,10);
stds1=zeros(10,10);
stds2=zeros(10,10);
stds3=zeros(10,10);
stds4=zeros(10,10);
stds5=zeros(10,10);
stds6=zeros(10,10);
stds7=zeros(10,10);
stds8=zeros(10,10);
avex=zeros(100,1);
stdx1=zeros(100,1);
stdx2=zeros(100,1);
stdx3=zeros(100,1);
stdx4=zeros(100,1);
stdx5=zeros(100,1);
stdx6=zeros(100,1);
stdx7=zeros(100,1);
stdx8=zeros(100,1);

bb4x = zeros(5,17);
bb5x = zeros(5,17);
bb6x = zeros(5,17);
bb7x = zeros(5,17);
bb8x = zeros(5,17);
stats = zeros(5,17);
%%name the files
pathx = 'f:\research\LiuWu\Database\DB';
pathfilex = ['01';'05';'10';'15';'17';'18';'20';'22';'27';'28';'29';'32';'33';'39';'42';'43';'45'];
pathfiles = strcat(pathx,pathfilex,'\datax1');
%% read data file and analyze
    % create discrete parameter values
    cos2 = zeros(10,1);
    cost = zeros(11,1);
    pre = zeros(11,1);
    cos2x = zeros(100,1);
    cos2x2 = zeros(100,1);
    costx = zeros(100,1);
    costx2 = zeros(100,1);
    coscost = zeros(100,1);
    cospre = zeros(100,1);
    prex = zeros(100,1);
    prex2 = zeros(100,1);
    
    % divide each parameter into 10 bins.
    for i = 1:10
        cos2(i) = (i-1)*0.1;
    end
    for i = 1:11
        cost(i) = (i-1)*0.1;
    end
    for i = 1:11
        %     pre(i) = (i-1)*0.20;
        pre(i) = (i-1)*0.15;
    end
%     ki = 1;
for ii = 1:size(pathfilex,1)
    fid = fopen(pathfiles(ii,:),'r');
    ki=1;
    tline = fgetl(fid); % after reading each line, move onto the next line.
    if ~ischar(tline), break, end
    while ~feof(fid) 
        tline = fgetl(fid); % after reading each line, move onto the next line.
        if ~ischar(tline), break, end
        [costwo pree1 pree2 aT bT cT xT yT zT mT nT] = strread(tline,'%f %f %f %f %f %f %f %f %f %f %f',1);
%         [v costwo pree preq costxx aT bT cT xT yT zT] = strread(tline,'%f %f %f %f %f %f %f %f %f %f %f',1);
        pree = sqrt(pree1*pree1+pree2*pree2);
%         pree = pree1*pree1+pree2*pree2;
        original(ki,:) = [sqrt(costwo) sqrt(pree) aT bT cT xT yT zT mT nT];
%         original(ki,:) = [sqrt(costwo) sqrt(costxx) sqrt(pree) aT bT cT xT yT zT];
        ki=ki+1;                  
    end % of while statement
    fclose(fid);
    clear fid;
    clear tline;

    % read data file done...
    %% group the data into 100 bins
    % classify the data according to the square of cosine theta(where theta is the angle
    % between the principal axis and treatment beam, and the cost function.
    % The cost function may exceed 1, in that case, the data points will not be
    % considered.
    inds1=ones(10,10);                   
    inds2=ones(10,10); 
    for k = 1:ki-1
        for i = 1:9
            for j = 1:10
%                 if original(k,1)>cos2(i) && original(k,1)<= cos2(i+1) && original(k,2)>cost(j) && original(k,2)<=cost(j+1) 
%                     sorted1(i,j,inds1(i,j)) = original(k,4);
%                     sorted2(i,j,inds1(i,j)) = original(k,5);
%                     sorted3(i,j,inds1(i,j)) = original(k,6);
%                     inds1(i,j)=inds1(i,j)+1;
%                 end
                if original(k,1)>cos2(i) && original(k,1)<= cos2(i+1) && original(k,2)>pre(j) && original(k,2)<=pre(j+1)
                    sorted4(i,j,inds2(i,j)) = original(k,6);
                    sorted5(i,j,inds2(i,j)) = original(k,7);
                    sorted6(i,j,inds2(i,j)) = original(k,8);
                    sorted7(i,j,inds2(i,j)) = original(k,9);
                    sorted8(i,j,inds2(i,j)) = original(k,10);
                    inds2(i,j)=inds2(i,j)+1;
                end
            end
        end
            % the cosine function has a limit.
        for j = 1:10
                %             if original(k,1)>cos2(10) && original(k,2)>cost(j) && original(k,2)<=cost(j+1) 
                %                 sorted1(10,j,inds1(10,j)) = original(k,4);
                %                 sorted2(10,j,inds1(10,j)) = original(k,5);
                %                 sorted3(10,j,inds1(10,j)) = original(k,6);                
                %                 inds1(10,j)=inds1(10,j)+1;
                %             end
            if original(k,1)>cos2(10) && original(k,2)>pre(j) && original(k,2)<=pre(j+1)             
                sorted4(10,j,inds2(10,j)) = original(k,6);
                sorted5(10,j,inds2(10,j)) = original(k,7);
                sorted6(10,j,inds2(10,j)) = original(k,8);
                sorted7(10,j,inds2(10,j)) = original(k,9);
                sorted8(10,j,inds2(10,j)) = original(k,10);
                inds2(10,j)=inds2(10,j)+1;
            end
        end      
    end
        for i = 1:10
            for j = 1:10        
                %             stds1(i,j)=std(sorted1(i,j,1:inds1(i,j)-1));
                %             stds2(i,j)=std(sorted2(i,j,1:inds1(i,j)-1));
                %             stds3(i,j)=std(sorted3(i,j,1:inds1(i,j)-1));
                stds4(i,j)=std(sorted4(i,j,1:inds2(i,j)-1));
                stds5(i,j)=std(sorted5(i,j,1:inds2(i,j)-1));
                stds6(i,j)=std(sorted6(i,j,1:inds2(i,j)-1));
                stds7(i,j)=std(sorted7(i,j,1:inds2(i,j)-1));
                stds8(i,j)=std(sorted8(i,j,1:inds2(i,j)-1));
                %             aves1(i,j)=mean(sorted1(i,j,1:inds1(i,j)-1));
                %             aves2(i,j)=mean(sorted2(i,j,1:inds1(i,j)-1));
                %             aves3(i,j)=mean(sorted3(i,j,1:inds1(i,j)-1));
                aves4(i,j)=mean(sorted4(i,j,1:inds2(i,j)-1));
                aves5(i,j)=mean(sorted5(i,j,1:inds2(i,j)-1));
                aves6(i,j)=mean(sorted6(i,j,1:inds2(i,j)-1));   
                aves7(i,j)=mean(sorted7(i,j,1:inds2(i,j)-1));
                aves8(i,j)=mean(sorted8(i,j,1:inds2(i,j)-1));            
            end
        end
        % flatten the stds and aves into one dimensional array, as well as the parameters into 1D array.
        for i = 1:100
            cos2x(i) = cos2(ceil(i/10))+0.05; % the square of cosine theta, where theta is the angle between the principal axis and treatment bean.   
            cos2x2(i) = cos2x(i)^2;
            if mod(i,10)>0 
                j = mod(i,10);
                costx(i) = cost(j)+0.05;        
                %         prex(i) = pre(j)+0.1;
                prex(i) = pre(j)+0.07;
                costx2(i)=costx(i)^2;
                prex2(i)=prex(i)^2;
                cospre(i) = cos2x(i)*prex(i);
                coscost(i) = cos2x(i)*costx(i);
            else
                j = 10;
                costx(i) = cost(10)+0.05;  
                %         prex(i) = pre(j)+0.1;
                prex(i) = pre(j)+0.07;
                costx2(i)=costx(i)^2;
                prex2(i)=prex(i)^2;
                cospre(i) = cos2x(i)*prex(i);
                coscost(i) = cos2x(i)*costx(i);
            end
            %     avex(i) = aves(ceil(i/10.0),j);
            %         stdx1(i) = stds1(ceil(i/10.0),j);            
            %         stdx2(i) = stds2(ceil(i/10.0),j); 
            %         stdx3(i) = stds3(ceil(i/10.0),j); 
            stdx4(i) = stds4(ceil(i/10.0),j);        
            stdx5(i) = stds5(ceil(i/10.0),j); 
            stdx6(i) = stds6(ceil(i/10.0),j);
            stdx7(i) = stds7(ceil(i/10.0),j);
            stdx8(i) = stds8(ceil(i/10.0),j);
        end
        % [X,Y]=meshgrid(cos2(1:10),cost(1:10));
        % figure(1)
        % surf(cos2(1:10),cost(1:10),stds1);
% colorbar;
% xlabel('|cos\theta)|', 'FontSize', 12, 'FontWeight', 'Bold');
% ylabel('sqrt(cost)', 'FontSize', 12, 'FontWeight', 'Bold');
% zlabel('standard deviation', 'FontSize', 12, 'FontWeight', 'Bold');
% % subplot(2,3,2);
% figure(2)
% surf(cos2(1:10),cost(1:10),stds2);
% colorbar;
% xlabel('|cos\theta)|', 'FontSize', 12, 'FontWeight', 'Bold');
% ylabel('sqrt(cost)', 'FontSize', 12, 'FontWeight', 'Bold');
% zlabel('standard deviation', 'FontSize', 12, 'FontWeight', 'Bold');
% figure(3);
% surf(cos2(1:10),cost(1:10),stds3);
% colorbar;
% xlabel('|cos\theta)|', 'FontSize', 12, 'FontWeight', 'Bold');
% ylabel('sqrt(cost)', 'FontSize', 12, 'FontWeight', 'Bold');
% zlabel('standard deviation', 'FontSize', 12, 'FontWeight', 'Bold');
% figure(4);
% surf(cos2(1:10),cost(1:10),stds4);
% colorbar;
% xlabel('|cos\theta)|', 'FontSize', 12, 'FontWeight', 'Bold');
% ylabel('sqrt(cost)', 'FontSize', 12, 'FontWeight', 'Bold');
% zlabel('standard deviation', 'FontSize', 12, 'FontWeight', 'Bold');
% figure(5);
% surf(cos2(1:10),cost(1:10),stds5);
% colorbar;
% xlabel('|cos\theta)|', 'FontSize', 12, 'FontWeight', 'Bold');
% ylabel('sqrt(cost)', 'FontSize', 12, 'FontWeight', 'Bold');
% zlabel('standard deviation', 'FontSize', 12, 'FontWeight', 'Bold');
% figure(6);
% surf(cos2(1:10),cost(1:10),stds6);
% colorbar;
% xlabel('|cos\theta)|', 'FontSize', 12, 'FontWeight', 'Bold');
% ylabel('sqrt(cost)', 'FontSize', 12, 'FontWeight', 'Bold');
% zlabel('standard deviation', 'FontSize', 12, 'FontWeight', 'Bold');
%     [bb1,bint1,r1,rint1,stats1] = regress(stdx1,[ones(100,1) cos2x cos2x2 costx] );
%     [bb2,bint2,r2,rint2,stats2] = regress(stdx2,[ones(100,1) cos2x cos2x2 costx] );
%     [bb3,bint3,r3,rint3,stats3] = regress(stdx3,[ones(100,1) cos2x cos2x2 costx] );
    % [bb3,bint3,r3,rint3,stats3] = regress(stdx3,[ones(100,1) cos2x cos2x2 costx coscost] );
    [bb4,bint4,r4,rint4,stats4] = regress(stdx4,[ones(100,1) cos2x cos2x2 prex prex2] );
    % [bb5,bint5,r5,rint5,stats5] = regress(stdx5,[ones(100,1) cos2x cos2x2 prex] );
    % [bb6,bint6,r6,rint6,stats6] = regress(stdx6,[ones(100,1) cos2x cos2x2 prex] );
    [bb5,bint5,r5,rint5,stats5] = regress(stdx5,[ones(100,1) cos2x cos2x2 prex prex2] );
    [bb6,bint6,r6,rint6,stats6] = regress(stdx6,[ones(100,1) cos2x cos2x2 prex prex2] );
    [bb7,bint7,r7,rint7,stats7] = regress(stdx7,[ones(100,1) cos2x cos2x2 prex prex2] );
    [bb8,bint8,r8,rint8,stats8] = regress(stdx8,[ones(100,1) cos2x cos2x2 prex prex2] );
% [bb5,bint5,r5,rint5,stats5] = regress(stdx5,[ones(100,1) cos2x cos2x2 prex cospre] );
% [bb6,bint6,r6,rint6,stats6] = regress(stdx6,[ones(100,1) cos2x cos2x2 prex cospre] );
    bb4x(:,ii) = bb4(:);
    bb5x(:,ii) = bb5(:);
    bb6x(:,ii) = bb6(:);
    bb7x(:,ii) = bb7(:);
    bb8x(:,ii) = bb8(:);
    stats(1,ii) = stats4(1);
    stats(2,ii) = stats5(1);
    stats(3,ii) = stats6(1);
    stats(4,ii) = stats7(1);
    stats(5,ii) = stats8(1);
end
pause(1);
