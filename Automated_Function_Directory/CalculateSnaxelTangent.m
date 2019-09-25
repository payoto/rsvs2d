function [varargout]=CalculateSnaxelTangent(varargin)
% include_SnakeParam
global CalculateSnaxelTangent_Handle
try
nOut=nargout(CalculateSnaxelTangent_Handle);
catch
include_SnakeParam
nOut=nargout(CalculateSnaxelTangent_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculateSnaxelTangent_Handle(varargin{:});
end
