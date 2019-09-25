function [varargout]=WritePolyStruct2Cell(varargin)
% include_PostProcessing
global WritePolyStruct2Cell_Handle
try
nOut=nargout(WritePolyStruct2Cell_Handle);
catch
include_PostProcessing
nOut=nargout(WritePolyStruct2Cell_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WritePolyStruct2Cell_Handle(varargin{:});
end
