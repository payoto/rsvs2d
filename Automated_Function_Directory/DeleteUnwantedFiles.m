function [varargout]=DeleteUnwantedFiles(varargin)
% include_Optimisation
global DeleteUnwantedFiles_Handle
nOut=nargout(DeleteUnwantedFiles_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DeleteUnwantedFiles_Handle(varargin{:});
end
