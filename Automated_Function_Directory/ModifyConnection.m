function [varargout]=ModifyConnection(varargin)
% include_NURBSEngine
global ModifyConnection_Handle
try
nOut=nargout(ModifyConnection_Handle);
catch
include_NURBSEngine
nOut=nargout(ModifyConnection_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ModifyConnection_Handle(varargin{:});
end
