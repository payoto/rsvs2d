function [varargout]=InverseDesign_ErrorTopo2(varargin)
global InverseDesign_ErrorTopo2_Handle
nOut=nargout(InverseDesign_ErrorTopo2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InverseDesign_ErrorTopo2_Handle(varargin{:});
end
