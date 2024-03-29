function popset(handle, prop)
% POPSET - pop values for a property from a value stack.
%
% POPSET(HANDLE, PROP) will restore a prevously HGPUSHED property value.
%

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

%  nargchk(2,2,'wrong number of arguments.');

if isempty(handle)
  return;
end

proplist = fieldnames(get(handle(1)));
prop = proplist{strcmpi(prop,proplist)};

appstr = [prop '_hgstack'];

for k=1:numel(handle)

  olds = getappdata(handle(k),appstr);

  if length(olds) <= 1
    warning(['Nothing left to pop for property ' prop '.']);
    continue;
  end

  set(handle(k),prop,olds{1});
  setappdata(handle(k),appstr,olds{2:end});

end

end
