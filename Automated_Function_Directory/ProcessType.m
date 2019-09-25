function [varargout]=ProcessType(varargin)
% include_Utilities
global ProcessType_Handle
try
nOut=nargout(ProcessType_Handle);
catch
include_Utilities
nOut=nargout(ProcessType_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProcessType_Handle(varargin{:});
end
