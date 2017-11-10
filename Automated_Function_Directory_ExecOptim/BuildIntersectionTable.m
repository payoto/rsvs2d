function [varargout]=BuildIntersectionTable(varargin)
global BuildIntersectionTable_Handle
nOut=nargout(BuildIntersectionTable_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildIntersectionTable_Handle(varargin{:});
end
