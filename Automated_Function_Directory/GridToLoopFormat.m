function [varargout]=GridToLoopFormat(varargin)
% include_Optimisation
global GridToLoopFormat_Handle
nOut=nargout(GridToLoopFormat_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GridToLoopFormat_Handle(varargin{:});
end
