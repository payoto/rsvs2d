function [varargout]=EditPLTHeader(varargin)
% OptimisationOutput
global EditPLTHeader_Handle
nOut=nargout(EditPLTHeader_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EditPLTHeader_Handle(varargin{:});
end
