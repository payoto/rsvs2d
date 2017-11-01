function [varargout]=WritePolyStruct2Cell(varargin)
% include_PostProcessing
global WritePolyStruct2Cell_Handle
nOut=nargout(WritePolyStruct2Cell_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WritePolyStruct2Cell_Handle(varargin{:});
end
