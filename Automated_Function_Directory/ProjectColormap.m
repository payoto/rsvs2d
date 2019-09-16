function [varargout]=ProjectColorMap(varargin)
% include_PostProcessing
global ProjectColorMap_Handle
nOut=nargout(ProjectColorMap_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProjectColorMap_Handle(varargin{:});
end
