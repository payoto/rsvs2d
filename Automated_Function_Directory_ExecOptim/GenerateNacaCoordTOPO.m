function [varargout]=GenerateNacaCoordTOPO(varargin)
global GenerateNacaCoordTOPO_Handle
nOut=nargout(GenerateNacaCoordTOPO_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateNacaCoordTOPO_Handle(varargin{:});
end
