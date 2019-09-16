function [varargout]=ContourNormal2(varargin)
% include_NURBSEngine
global ContourNormal2_Handle
nOut=nargout(ContourNormal2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ContourNormal2_Handle(varargin{:});
end
