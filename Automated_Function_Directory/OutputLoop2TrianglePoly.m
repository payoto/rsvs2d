function [varargout]=OutputLoop2TrianglePoly(varargin)
% include_PostProcessing
global OutputLoop2TrianglePoly_Handle
nOut=nargout(OutputLoop2TrianglePoly_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OutputLoop2TrianglePoly_Handle(varargin{:});
end
