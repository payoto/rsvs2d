function [varargout]=REFINE_contcurvevolnoedge(varargin)
global REFINE_contcurvevolnoedge_Handle
nOut=nargout(REFINE_contcurvevolnoedge_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_contcurvevolnoedge_Handle(varargin{:});
end
