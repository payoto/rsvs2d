function [varargout]=InverseDesign(varargin)
global InverseDesign_Handle
nOut=nargout(InverseDesign_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InverseDesign_Handle(varargin{:});
end
