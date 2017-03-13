function [varargout]=REFINE_contcurve(varargin)
global REFINE_contcurve_Handle
nOut=nargout(REFINE_contcurve_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_contcurve_Handle(varargin{:});
end
