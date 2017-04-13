function [varargout]=DisplacementsOutput(varargin)
global DisplacementsOutput_Handle
nOut=nargout(DisplacementsOutput_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DisplacementsOutput_Handle(varargin{:});
end
