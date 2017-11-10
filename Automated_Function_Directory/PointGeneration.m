function [varargout]=PointGeneration(varargin)
% include_SnakeParam
global PointGeneration_Handle
nOut=nargout(PointGeneration_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PointGeneration_Handle(varargin{:});
end
