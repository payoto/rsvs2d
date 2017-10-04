function [varargout]=PlotGradients(varargin)
% OptimisationOutput
global PlotGradients_Handle
nOut=nargout(PlotGradients_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotGradients_Handle(varargin{:});
end
