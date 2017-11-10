function [varargout]=PlotGradientsRefine(varargin)
% OptimisationOutput
global PlotGradientsRefine_Handle
nOut=nargout(PlotGradientsRefine_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=PlotGradientsRefine_Handle(varargin{:});
end
