function [snaxmove]=BuildMovementStructures(snaxel,snakposition,snaxmode,sensSnax,volumefraction)
  
    
    nSnax=length(snaxel);
    
    [snaxmove]=BuildSnaxMoveTemplate(nSnax,nSnax,nSnax);
    
    
end


function [snaxmove]=BuildSnaxMoveTemplate(nSnax,nEdge,nVertex)
    
    snaxmove=struct('snax',struct([]),'vertex',struct([]),'edge',struct([]),...
        'support',struct('nSnax',[],'maxL',[]));
    
    snaxTemp=struct('index',[],'rN',[],'rT',[],'d',[],'sens',[],'posL',[]);
    
    edgeTemp=struct('posL',[],'coord',[],'normal',[]);
    vertTemp=edgeTemp;
    
    snaxmove.snax=repmat(snaxTemp,[1,nSnax]);
    snaxmove.vertex=repmat(edgeTemp,[1,nVertex]);
    snaxmove.edge=repmat(vertTemp,[1,nEdge]);
    
end


function [snaxmove]=CalculateVertexData(snaxmove,snaxel,snakposition,sensSnax)
    
    
    
    for ii=1:length(snakposition)
       snaxmove.snax(ii).index=snakposition(ii).index;
       snaxmove.snax(ii).d=sqrt(sum(snakposition(ii).vectornotnorm.^2));
       snaxmove.snax(ii).sens=sensSnax(ii,:);
       
       snaxmove.snax(ii).index=snakposition(ii).index;
       snaxmove.snax(ii).index=snakposition(ii).index; 
        
    end
    
    
    
    
end