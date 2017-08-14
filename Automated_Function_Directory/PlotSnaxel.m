function [varargout]=PlotSnaxel(varargin)
% include_CheckResultsLight
global PlotSnaxel_Handle
nOut=nargout(PlotSnaxel_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotSnaxel_Handle(varargin{:});
end
