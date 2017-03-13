function [varargout]=Refinement_edgecross(varargin)
global Refinement_edgecross_Handle
nOut=nargout(Refinement_edgecross_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=Refinement_edgecross_Handle(varargin{:});
end
