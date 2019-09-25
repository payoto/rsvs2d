function [varargout]=ModifUnstructured(varargin)
% include_SnakeParam
global ModifUnstructured_Handle
try
nOut=nargout(ModifUnstructured_Handle);
catch
include_SnakeParam
nOut=nargout(ModifUnstructured_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ModifUnstructured_Handle(varargin{:});
end
