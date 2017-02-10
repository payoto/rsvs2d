function [varargout]=REFINE_contour(varargin)
global REFINE_contour_Handle
nOut=nargout(REFINE_contour_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=REFINE_contour_Handle(varargin{:});
end
