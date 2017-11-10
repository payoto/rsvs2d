function [varargout]=GenerateSurrogateStructure(varargin)
global GenerateSurrogateStructure_Handle
nOut=nargout(GenerateSurrogateStructure_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateSurrogateStructure_Handle(varargin{:});
end
