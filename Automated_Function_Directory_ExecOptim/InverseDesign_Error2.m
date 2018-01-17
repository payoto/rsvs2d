function [varargout]=InverseDesign_Error2(varargin)
% InverseDesign_Error
global InverseDesign_Error2_Handle
nOut=nargout(InverseDesign_Error2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InverseDesign_Error2_Handle(varargin{:});
end
