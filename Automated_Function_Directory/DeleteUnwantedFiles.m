function [varargout]=DeleteUnwantedFiles(varargin)
% include_Optimisation
global DeleteUnwantedFiles_Handle
try
nOut=nargout(DeleteUnwantedFiles_Handle);
catch
include_Optimisation
nOut=nargout(DeleteUnwantedFiles_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DeleteUnwantedFiles_Handle(varargin{:});
end
