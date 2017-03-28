function [varargout]=REFINE_contcurvevol(varargin)
global REFINE_contcurvevol_Handle
nOut=nargout(REFINE_contcurvevol_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_contcurvevol_Handle(varargin{:});
end
