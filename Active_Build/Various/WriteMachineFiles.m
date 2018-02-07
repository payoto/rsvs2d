function [machinedat]=WriteMachineFiles(nWorker,writeDir)
    
    
    %     load('jobdat.mat','pbsdat')
    %     pbsdat.numppnnode=16;
    %
    %     pbsdat.numppnnode=ones([1,pbsdat.numnode])*(pbsdat.numppnnode);
    %     pbsdat.numppnnode(1)=pbsdat.numppnnode(1)-1;
    %     pbsdat.numppntotal=sum(pbsdat.numppnnode);
    %     pbsdat.jobid='6567677-1.master.cm.cluster';
    
    comStr=computer;
    if strcmp(comStr(1:2),'PC')
        pbsdat.jobid=int2str(randi(9999));
        [~,pbsdat.nodelist]=system('hostname');
        pbsdat.nodelist=regexprep(pbsdat.nodelist,'\n','');
        pbsdat.numnode=1; % Number of unique nodes
        [~,numCores]=system('WMIC CPU Get NumberOfCores');
        numCores=str2double(regexprep(numCores,'NumberOfCores',''));
        matInd=1;
        pbsdat.numppnnode=ones([1,pbsdat.numnode])*numCores; % Number of processor per node
        pbsdat.numppnnode(matInd)=pbsdat.numppnnode(matInd)-1;
        pbsdat.numppntotal=sum(pbsdat.numppnnode); % total number of processors
    else
        [pbsdat]=CollectPBSData();
    end
    
    [machinedat]=BuildMachineDat(nWorker,pbsdat);
    marker=regexp(pbsdat.jobid,'^.*[0-9]','match');
    marker=marker{1};
    [machinedat]=GenerateMachineFiles(writeDir,machinedat,marker);
end

function [machinedat]=GenerateMachineFiles(writeDir,machinedat,marker)
    
    
    for ii=1:numel(machinedat)
        machinedat(ii).file=[writeDir,filesep,'mpihosts_',int2str(ii),'_',marker];
        fid=fopen(machinedat(ii).file,'w');
        if fid==0
            error('Generate Machine Files failed')
        end
        for jj=1:size(machinedat(ii).node,1)
            
            fprintf(fid,'%s slots=%i max_slots=%i\n',machinedat(ii).node(jj,:),...
                machinedat(ii).nodeslots(jj,:),machinedat(ii).nodeslots(jj,:));
        end
        fclose(fid);
    end
    
end

function [machinedat]=BuildMachineDat(nWorker,pbsdat)
    
    
    machinedat=repmat(struct('node','','nodeslots',[],'slots',[],'file',''),[1,nWorker]);
    
    
    meanPPNNode=mean(pbsdat.numppnnode);
    iWork=0;
    actNumNode=pbsdat.numppnnode;
    for ii=1:pbsdat.numnode
        meanWorkerNode=(nWorker-iWork)/(pbsdat.numnode-(ii-1));
        if meanPPNNode>pbsdat.numppnnode(ii)
            targWorkerNode=floor(meanWorkerNode);
        else
            targWorkerNode=ceil(meanWorkerNode);
        end
        while targWorkerNode>0
            iWork=iWork+1;
            machinedat(iWork).node=pbsdat.nodelist(ii,:);
            machinedat(iWork).nodeslots=pbsdat.numppnnode(ii);
            ppnLeft=actNumNode(ii)/targWorkerNode;
            
            machinedat(iWork).slots=round(ppnLeft);
            actNumNode(ii)=actNumNode(ii)-machinedat(iWork).slots;
            targWorkerNode=targWorkerNode-1;
        end
    end
end

function [pbsdat]=CollectPBSData()
    % Collects PBS Data
    
    [~,jobid]=system('echo $PBS_JOBID');
    [~,jobhost]=system('echo $PBS_O_HOST');
    
    [~,jobNodeFile]=system('echo $PBS_NODEFILE');
    [~,nodeListFull]=system('cat $PBS_NODEFILE');
    [~,nodeList]=system('cat $PBS_NODEFILE | uniq');
    
    [~,numNodes]=system('echo $PBS_NUM_NODES');
    [~,numPPNNode]=system('echo $PBS_NUM_PPN');
    
    matNode=regexp(nodeListFull,'\n','split');
    
    nodeListCell=regexp(nodeList,'\n','split');
    nodeListCell=nodeListCell(~cellfun(@isempty,nodeListCell));
    
    pbsdat.matnode=matNode{1}; % Node Matlab is running on
    
    matInd=find(~cellfun(@isempty,regexp(nodeListCell,pbsdat.matnode)));
    
    pbsdat.jobid=jobid;
    pbsdat.nodelist=char(nodeListCell{:}); % List of unique nodes
    pbsdat.numnode=str2double(numNodes); % Number of unique nodes
    pbsdat.numppnnode=ones([1,pbsdat.numnode])*str2double(numPPNNode); % Number of processor per node
    pbsdat.numppnnode(matInd)=pbsdat.numppnnode(matInd)-1;
    pbsdat.numppntotal=sum(pbsdat.numppnnode); % total number of processors
    %save('jobdat.mat')
end