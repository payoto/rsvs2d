function [varargout]=GenerateNURBSArea(varargin)
% include_NURBSEngine
global GenerateNURBSArea_Handle
try
nOut=nargout(GenerateNURBSArea_Handle);
catch
include_NURBSEngine
nOut=nargout(GenerateNURBSArea_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateNURBSArea_Handle(varargin{:});
end
