function [varargout]=ClosedLoopProperties(varargin)
% include_Utilities
global ClosedLoopProperties_Handle
nOut=nargout(ClosedLoopProperties_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ClosedLoopProperties_Handle(varargin{:});
end
