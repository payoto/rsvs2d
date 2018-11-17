function [machinedat]=WriteMachineFiles(nWorker,writeDir,node1ReserveNum)
    
    
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
        pbsdat.numppnnode(matInd)=pbsdat.numppnnode(matInd)-node1ReserveNum;
        pbsdat.numppntotal=sum(pbsdat.numppnnode); % total number of processors
    else
        [~,hostName]=system('whichbluecrystal');
        if ~isempty(regexp(hostName,'4', 'once'))
            [pbsdat]=CollectSLURMData(node1ReserveNum);
        elseif ~isempty(regexp(hostName,'3', 'once'))
            [pbsdat]=CollectPBSData(node1ReserveNum);
        end
        
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
            
            fprintf(fid,'%s slots=%i max-slots=%i\n',machinedat(ii).node(jj,:),...
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
            slotindices = actNumNode(ii)-round(ppnLeft)+1:actNumNode(ii);
            cpuset = strjoin(strsplit(num2str(slotindices),' '),',');
            machinedat(iWork).cpuset=cpuset;
            actNumNode(ii)=actNumNode(ii)-machinedat(iWork).slots;
            targWorkerNode=targWorkerNode-1;
        end
    end
end

function [pbsdat]=CollectPBSData(node1ReserveNum)
    % Collects PBS Data
    %  node1ReserveNum - no. of slots to hold back on main node (for MATLAB etc.)
    
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
    pbsdat.numppnnode(matInd)=pbsdat.numppnnode(matInd)-node1ReserveNum;
    pbsdat.numppntotal=sum(pbsdat.numppnnode); % total number of processors
    %save('jobdat.mat')
end

function [pbsdat]=CollectSLURMData(node1ReserveNum)
    % Collects PBS Data
    %  node1ReserveNum - no. of slots to hold back on main node (for MATLAB etc.)
    
    [~,jobid]=system('echo $SLURM_JOBID');
    
    [~,nodeListFull]=system('echo $SLURM_JOB_NODELIST'); % format is compute[092,387]
    
    
    [~,numNodes]=system('echo $SLURM_NNODES');
    [~,numPPNNode]=system('echo $SLURM_NTASKS_PER_NODE');
    
    
    nodeList=SlurmToNodeList(nodeListFull);
    [~,host]=system('hostname');
    matNode=regexp(host,'\.','split');
    for ii=1:size(nodeList,1)
        nodeListCell{ii}=nodeList(ii,:);
    end
    nodeListCell=nodeListCell(~cellfun(@isempty,nodeListCell));
    
    pbsdat.matnode=matNode{1}; % Node Matlab is running on
    
    matInd=find(~cellfun(@isempty,regexp(nodeListCell,pbsdat.matnode)));
    
    pbsdat.jobid=jobid;
    pbsdat.nodelist=nodeList;
    pbsdat.numnode=str2double(numNodes); % Number of unique nodes
    pbsdat.numppnnode=ones([1,pbsdat.numnode])*str2double(numPPNNode); % Number of processor per node
    pbsdat.numppnnode(matInd)=pbsdat.numppnnode(matInd)-node1ReserveNum;
    pbsdat.numppntotal=sum(pbsdat.numppnnode); % total number of processors
    %save('jobdat.mat')
end

function [nodelist]=SlurmToNodeList(slurmOut)
    % takes the output of a echo $SLURM_JOB_NODELIST
    % and turns it into a list of nodes as a character array
    
    slurmCells=regexp(slurmOut,'\]','split');
    slurmCells=regexprep(slurmCells,'^,','');
    nodeCells=cell(0);
    kk=0;
    for ii=1:numel(slurmCells)
        if ~isempty(regexp(slurmCells{ii},'\S','once'))
            sepPoint=regexp(slurmCells{ii},'\[.*');
            if isempty(sepPoint)
                nodeCells{kk+1}=slurmCells{ii};
                kk=kk+1;
            else
                baseName=slurmCells{ii}(1:sepPoint-1);
                
                numArray=regexp(slurmCells{ii}(sepPoint+1:end),',','split');
                
                for jj=1:numel(numArray)
                    if numel(numArray{jj})==3
                        nodeCells{kk+1}=[baseName,numArray{jj}];
                        kk=kk+1;
                    else % case where the num appears as ###-###
                        intArray=cellfun(@str2double,regexp(numArray{jj},'-','split'));
                        for ll=intArray(1):intArray(2)
                            nodeCells{kk+1}=[baseName,num2str(ll,'%3i')];
                            kk=kk+1;
                        end
                    end
                end
            end
        end
    end
    nodeCells=regexprep(nodeCells,'\s','');
    nodeCells=deblank(nodeCells(~cellfun(@isempty,nodeCells)));
    nodelist=char(nodeCells{:});
    
end



%{
PARTITION     AVAIL  TIMELIMIT  NODES  STATE NODELIST
veryshort        up    6:00:00      5 drain* compute[223,249,259,271,485]
veryshort        up    6:00:00     21   comp compute[161,167,179,186,193,205,227,268,299,337,342,370,386,397,401,406,412-413,425,432,460]
veryshort        up    6:00:00     14  drain compute[395,418,433-444]
veryshort        up    6:00:00      3   resv compute[082,123,185]
veryshort        up    6:00:00     20    mix compute[119-120,122,125,158,170,176,188,191,239,243-244,246-247,252,305-306,341,492,522]
veryshort        up    6:00:00    388  alloc compute[076-081,083-097,099-118,121,124,126-157,159-160,162-166,168-169,171-175,177-178,180-184,187,189-190,192,194-204,206-222,224-226,228-238,240-242,245,248,250-251,253-258,260-267,269-270,272-298,300-304,307,309-336,338-340,343-369,371-385,387-394,396,398-400,402-405,407-411,414-417,419-424,426-430,445-459,461-484,486-491,493-521,523-525],highmem[10-13]
veryshort        up    6:00:00      3   idle compute[098,308,431]
test             up    1:00:00      5 drain* compute[223,249,259,271,485]
test             up    1:00:00     21   comp compute[161,167,179,186,193,205,227,268,299,337,342,370,386,397,401,406,412-413,425,432,460]
test             up    1:00:00     14  drain compute[395,418,433-444]
test             up    1:00:00      3   resv compute[082,123,185]
test             up    1:00:00     20    mix compute[119-120,122,125,158,170,176,188,191,239,243-244,246-247,252,305-306,341,492,522]
test             up    1:00:00    388  alloc compute[076-081,083-097,099-118,121,124,126-157,159-160,162-166,168-169,171-175,177-178,180-184,187,189-190,192,194-204,206-222,224-226,228-238,240-242,245,248,250-251,253-258,260-267,269-270,272-298,300-304,307,309-336,338-340,343-369,371-385,387-394,396,398-400,402-405,407-411,414-417,419-424,426-430,445-459,461-484,486-491,493-521,523-525],highmem[10-13]
test             up    1:00:00      3   idle compute[098,308,431]
cpu_test         up    1:00:00      1    mix compute075
cpu_test         up    1:00:00      1   idle compute074
cpu*             up 14-00:00:0      5 drain* compute[223,249,259,271,485]
cpu*             up 14-00:00:0     21   comp compute[161,167,179,186,193,205,227,268,299,337,342,370,386,397,401,406,412-413,425,432,460]
cpu*             up 14-00:00:0     14  drain compute[395,418,433-444]
cpu*             up 14-00:00:0      3   resv compute[082,123,185]
cpu*             up 14-00:00:0     20    mix compute[119-120,122,125,158,170,176,188,191,239,243-244,246-247,252,305-306,341,492,522]
cpu*             up 14-00:00:0    384  alloc compute[076-081,083-097,099-118,121,124,126-157,159-160,162-166,168-169,171-175,177-178,180-184,187,189-190,192,194-204,206-222,224-226,228-238,240-242,245,248,250-251,253-258,260-267,269-270,272-298,300-304,307,309-336,338-340,343-369,371-385,387-394,396,398-400,402-405,407-411,414-417,419-424,426-430,445-459,461-484,486-491,493-521,523-525]
cpu*             up 14-00:00:0      3   idle compute[098,308,431]
hmem             up 14-00:00:0      5  alloc highmem[10-14]
hmem             up 14-00:00:0      3   idle highmem[15-17]
gpu              up 7-00:00:00      1  down* gpu14
gpu              up 7-00:00:00      2   comp gpu[09,24]
gpu              up 7-00:00:00     21    mix gpu[01-08,10-13,15-23]
gpu              up 7-00:00:00      3   idle gpu[25-27]
gpu_veryshort    up    1:00:00      2    mix gpu[28-29]
gpu_veryshort    up    1:00:00      3   idle gpu[30-32]
serial           up 3-00:00:00      1   comp compute070
serial           up 3-00:00:00      4    mix compute[068-069,071-072]
serial           up 3-00:00:00      1   idle compute073
dcv              up 14-00:00:0      1    mix bc4vis1
%}
