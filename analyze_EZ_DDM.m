% cl968@georgetown.edu

% define some variables
loads = {'OneBack','TwoBack','ThreeBack','FourBack'};
blocks = {'A','B','C','D','E'};
targets = {'NonHub','Hub'};

% define 
% number 
% of subjects 
nSubjects = 24;

% sweep
% through
% subjects
for s = 1:nSubjects
    
    % sweep
    % through
    % targets
    for t = 1:length(targets)
        
        % read in Excel data
        warning off; % suppress annoying warning
        data = dataset('XLSFile',['eprime/Nback-S'...
            num2str(s) '-' targets{t} '.xlsx']);
        
        % populate CRESP field
        % with correct answers
        for i = 1:length(data)
            if strcmp(data.Condition(i),'FILL')
                data.stimulus_CRESP(i)= 1;
            else
                data.stimulus_CRESP(i)= 2;
            end
        end
        
        % preallocate 
        % temporary variable 
        temp = zeros(length(data),3);
        
        % Note: rows in variable temp are individual trials. Column 1 is 
        % load conditions, column 2 is acc. (1=correct, 0=incorrect),
        % and column 3 holds response times in ms 

        % identify
        % incorrect
        % responses
        for i = 1:length(data)
            
            if data.stimulus_RESP(i)~=data.stimulus_CRESP(i)
                temp(i,2)=0;
            else
                temp(i,2)=1;
            end
            
        end
        
        % preallocate
        load_idx = zeros(length(data),...
            length(loads));
        
        % construct load index
        for i = 1:length(loads)
            
            % sweep blocks
            for ii = 1:length(blocks)
                
                % concatenate
                trials = data.([loads{i} blocks{ii}]);
                
                % sweep individual trials
                for iii = 1:length(trials)
                    
                    if ~strcmp(trials{iii,1},'NULL')
                        load_idx(iii,i)=i;
                    end
                    
                end
                
            end
            
        end
                     
        % conditions
        temp(:,1) = sum(load_idx,2); 
        
        % load in response times
        temp(:,3) = data.stimulus_RT;
        
        % remove fast responses 
        temp(temp(:,3)<=300,:)=[];
        
        % input parameters
        rt = temp(:,3); % response times
        acc = temp(:,2); % accuracy
        
        % sweep
        % through
        % loads
        for i = 1:length(loads)
            
            % load index
            lidx = find(temp(:,1)==i);
            
            % log input parameters
            ddm.(['l' num2str(i)]).acc(s,t) = mean(acc(lidx));
            ddm.(['l' num2str(i)]).vrt(s,t) = var(rt(temp(lidx,2)==1));
            ddm.(['l' num2str(i)]).mrt(s,t) = mean(rt(temp(lidx,2)==1));
            
            % log output parameters 
            [ddm.(['l' num2str(i)]).v(s,t),ddm.(['l' num2str(i)]).a(s,t),ddm.(['l' num2str(i)]).ter(s,t)] = ...
                calc_ddm(mean(acc(lidx)),var(rt(temp(lidx,2)==1)),mean(rt(temp(lidx,2)==1)));
            
        end
           
    end
    
end

% specify
% order of 
% diffusion
% parameters
fn = {'v','a','ter',...
    'acc','mrt','vrt'};

% preallocate 
corr_mat = ...
    zeros(6,6);

% sweep 
% through
% parameters 
for i = 1:length(fn)
    
    a=[]; % temp
    
    % sweep 
    % through loads 
    for w = 1:4
        a = [a;ddm.(['l' num2str(w)]).(fn{i})(:,1);ddm.(['l' num2str(w)]).(fn{i})(:,2)];
    end
    
    % sweep 
    % through parameters 
    for ii = 1:length(fn)
        
        b=[]; % temp
        
        % sweep
        % through loads 
        for w = 1:4
            b = [b;ddm.(['l' num2str(w)]).(fn{ii})(:,1);ddm.(['l' num2str(w)]).(fn{ii})(:,2)];
        end
        
        % log corr
        if i~=ii
            corr_mat(i,ii)=corr(a,b);
        end
        
    end
    
end

% remove upper half of matrix 
corr_mat(triu(ones(6),1) > 0)=0;

% create colormap 
bwr = b2r(min(min(corr_mat)),...
    max(max(corr_mat)));

% plot matrix
subplot(1,2,1)
imagesc(corr_mat);
colormap(bwr); hold
set(gca,'TickLength',[0 0],'XTickLabel',fn); % set x-axis
set(gca,'TickLength',[0 0],'YTickLabel',fn); % set x-axis
box('off')

% add 
% grid
hold on;
for i = 1:(11+1)    
     plot([.5,(i+.5)],[i-.5,i-.5],'k-');
     plot([i-.5,i-.5],[(-1.5+i),12.5],'k-'); 
end

% preallocate
N = zeros(nSubjects,length(loads)); % nonhubs
H = zeros(nSubjects,length(loads)); % hubs

% sweep through
% participants
for i = 1:nSubjects
    for ii = 1:length(loads)
        N(i,ii)=ddm.(['l' num2str(ii)]).(fn{1})(i,1);
        H(i,ii)=ddm.(['l' num2str(ii)]).(fn{1})(i,2);
    end
end

% preallocate
M = [];
SE = [];

% create
% plot variable
for i = 1:length(loads)
    M(i,1)=mean(N(:,i));
    M(i,2)=mean(H(:,i));
    SE(i,1)=std(N(:,i))/sqrt(nSubjects);
    SE(i,2)=std(H(:,i))/sqrt(nSubjects);
end

% plot 
% v bar
% graphs
subplot(1,2,2)
bar(M); hold;
errorbar([.85 1.15 1.85 2.15 2.85 3.15 3.85 4.15],...
    [M(1,:) M(2,:) M(3,:) M(4,:)],...
    [SE(1,:) SE(2,:) SE(3,:) SE(4,:)],'.','Color','k');
R = [H(:,1);N(:,1); H(:,2); N(:,2); H(:,3); N(:,3); H(:,4); N(:,4)];
set(gca,'TickLength',[0 0]); % set x-axis
ylabel(fn{1},'FontSize',14);
xlabel('Cognitive Load','FontSize',14);


% sweep
% through
% other
% parameters
for p = 2:3
    
    % preallocate
    N = zeros(nSubjects,length(loads)); % nonhubs
    H = zeros(nSubjects,length(loads)); % hubs
    
    % sweep through
    % participants
    for i = 1:nSubjects
        for ii = 1:length(loads)
            N(i,ii)=ddm.(['l' num2str(ii)]).(fn{p})(i,1);
            H(i,ii)=ddm.(['l' num2str(ii)]).(fn{p})(i,2);
        end
    end
    
    % preallocate
    M = [];
    SE = [];
    
    % create
    % plot variable
    for i = 1:length(loads)
        M(i,1)=mean(N(:,i));
        M(i,2)=mean(H(:,i));
        SE(i,1)=std(N(:,i))/sqrt(nSubjects);
        SE(i,2)=std(H(:,i))/sqrt(nSubjects);
    end
    
    % plot bar
    % graphs
    subplot(1,2,p-1)
    bar(M); hold;
    errorbar([.85 1.15 1.85 2.15 2.85 3.15 3.85 4.15],...
        [M(1,:) M(2,:) M(3,:) M(4,:)],...
        [SE(1,:) SE(2,:) SE(3,:) SE(4,:)],'.','Color','k');
    R = [H(:,1);N(:,1); H(:,2); N(:,2); H(:,3); N(:,3); H(:,4); N(:,4)];
    set(gca,'TickLength',[0 0]); % set x-axis
    ylabel(fn{p},'FontSize',14);
    xlabel('Cognitive Load','FontSize',14);
    
end





