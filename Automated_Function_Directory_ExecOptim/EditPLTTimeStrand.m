function [varargout]=EditPLTTimeStrand(varargin)
% OptimisationOutput
global EditPLTTimeStrand_Handle
nOut=nargout(EditPLTTimeStrand_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EditPLTTimeStrand_Handle(varargin{:});
end
