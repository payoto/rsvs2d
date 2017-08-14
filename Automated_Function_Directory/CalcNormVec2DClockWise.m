function [varargout]=CalcNormVec2DClockWise(varargin)
% include_EdgeInformation
global CalcNormVec2DClockWise_Handle
nOut=nargout(CalcNormVec2DClockWise_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalcNormVec2DClockWise_Handle(varargin{:});
end
