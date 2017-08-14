function [varargout]=WriteFullInfoProfile(varargin)
% include_PostProcessing
global WriteFullInfoProfile_Handle
nOut=nargout(WriteFullInfoProfile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WriteFullInfoProfile_Handle(varargin{:});
end
