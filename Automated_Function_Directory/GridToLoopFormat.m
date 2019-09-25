function [varargout]=GridToLoopFormat(varargin)
% include_Optimisation
global GridToLoopFormat_Handle
try
nOut=nargout(GridToLoopFormat_Handle);
catch
include_Optimisation
nOut=nargout(GridToLoopFormat_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GridToLoopFormat_Handle(varargin{:});
end
