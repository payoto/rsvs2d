function [varargout]=REFINE_contcurvescale(varargin)
global REFINE_contcurvescale_Handle
nOut=nargout(REFINE_contcurvescale_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_contcurvescale_Handle(varargin{:});
end
