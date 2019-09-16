function [varargout]=ModifyConnection(varargin)
% include_NURBSEngine
global ModifyConnection_Handle
nOut=nargout(ModifyConnection_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ModifyConnection_Handle(varargin{:});
end
