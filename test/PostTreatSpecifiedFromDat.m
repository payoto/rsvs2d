function []=PostTreatSpecifiedFromDat(fileName,nFiles)

InitialiseSnakeFlow

for ii=1:nFiles
	fid=fopen(fileName,'r')
	for jj=1:ii
		str=fgetl(fid);
	end
	load(str)
	kk=1;
	while ~isempty(optimstruct(kk).population(1).objective)
		kk=kk+1;
	end
	optimstruct(kk:end)=[];
	posPath=regexp(str,'/');
	
	PostTreatIncomplete(str(1:(posPath(end)-1)),[],optimstruct)
end

end
