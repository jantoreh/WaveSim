function Snew=mkdspec(S,D,plotflag)

% Make a directional spectrum combining frequency and directionally
% dependent S and D.
%
% Plotflag is still unused
  
if isempty(S)
  error('Input spectrum is empty.')
end
if nargin<2 || isempty(D)
  error('Input spreading is empty.')
end
if nargin<3
  plotflag=0;
end

Snew = struct();


if isfield(D,'note')
  Snew.note=[Snew.note,'; ',D.note];
else
  Snew.note=[];
end

Snew.date  = datestr(now);
Snew.theta = D.theta;

if ~isfield(D,'w') || isempty(D.w) || length(D.w)==1 % Not frequency dependent
  Snew.S = D.S(:)*S.S(:)';
else
  if length(S.w)~=length(D.w)
    error('Frequency in S and D must be identical.')
  elseif any(abs(S.w-D.w) >= 1e-10)
    error('Frequency in S and D must be identical.')
  end    
  Snew.S=D.S.*S.S(:,ones(1,length(D.theta)))';
end

if plotflag == 1
   % Plot spectrum 
end
