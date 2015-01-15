function [states] = stateLabels(nStates,len)
%
% Donada una longitud de cadena 'len' i el nombre d'estats 'nStates',
% genera la cadena 'sates' amb les etiquetes d'estat de cada frame,
% dividint la sequencia a parts iguals.

states=ones(1,len);
gap=round(len/nStates);

for i=2:nStates
    ini=gap*(i-1)+1;
    states(1,ini:ini+gap)=i;
end
states=states(1:len);

end

