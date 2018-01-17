function [varargout]=GenerateNacaCoord(varargin)
% InverseDesign_Error
global GenerateNacaCoord_Handle
nOut=nargout(GenerateNacaCoord_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateNacaCoord_Handle(varargin{:});
end
