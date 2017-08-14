function [varargout]=CalculateSnaxelTangent(varargin)
% include_SnakeParam
global CalculateSnaxelTangent_Handle
nOut=nargout(CalculateSnaxelTangent_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculateSnaxelTangent_Handle(varargin{:});
end
