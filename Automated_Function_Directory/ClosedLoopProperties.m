function [varargout]=ClosedLoopProperties(varargin)
% include_Utilities
global ClosedLoopProperties_Handle
try
nOut=nargout(ClosedLoopProperties_Handle);
catch
include_Utilities
nOut=nargout(ClosedLoopProperties_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ClosedLoopProperties_Handle(varargin{:});
end
