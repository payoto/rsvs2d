function [varargout]=GeneratePopulationStruct(varargin)
global GeneratePopulationStruct_Handle
nOut=nargout(GeneratePopulationStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GeneratePopulationStruct_Handle(varargin{:});
end
