function [varargout]=TecplotPortion_Init(varargin)
% OptimisationOutput
global TecplotPortion_Init_Handle
nOut=nargout(TecplotPortion_Init_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TecplotPortion_Init_Handle(varargin{:});
end
