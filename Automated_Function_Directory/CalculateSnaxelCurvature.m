function [varargout]=CalculateSnaxelCurvature(varargin)
global CalculateSnaxelCurvature_Handle
nOut=nargout(CalculateSnaxelCurvature_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculateSnaxelCurvature_Handle(varargin{:});
end
