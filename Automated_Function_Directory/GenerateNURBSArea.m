function [varargout]=GenerateNURBSArea(varargin)
% include_NURBSEngine
global GenerateNURBSArea_Handle
nOut=nargout(GenerateNURBSArea_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateNURBSArea_Handle(varargin{:});
end
