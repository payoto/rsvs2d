function [varargout]=GenerateNACACoordOrient(varargin)
global GenerateNACACoordOrient_Handle
nOut=nargout(GenerateNACACoordOrient_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateNACACoordOrient_Handle(varargin{:});
end
