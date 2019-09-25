function [varargout]=FindDir(varargin)
% include_Optimisation
global FindDir_Handle
try
nOut=nargout(FindDir_Handle);
catch
include_Optimisation
nOut=nargout(FindDir_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindDir_Handle(varargin{:});
end
