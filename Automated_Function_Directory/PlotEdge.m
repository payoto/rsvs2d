function [varargout]=PlotEdge(varargin)
% include_CheckResultsLight
global PlotEdge_Handle
nOut=nargout(PlotEdge_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotEdge_Handle(varargin{:});
end
