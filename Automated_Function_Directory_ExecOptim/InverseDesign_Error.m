function [varargout]=InverseDesign_Error(varargin)
% include_InverseDesign_Error
global InverseDesign_Error_Handle
nOut=nargout(InverseDesign_Error_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InverseDesign_Error_Handle(varargin{:});
end
