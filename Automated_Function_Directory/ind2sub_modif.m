function [varargout]=ind2sub_modif(varargin)
global ind2sub_modif_Handle
nOut=nargout(ind2sub_modif_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ind2sub_modif_Handle(varargin{:});
end
