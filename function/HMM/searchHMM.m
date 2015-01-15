function [ seqLikelihoods ] = searchHMM(lenMin, lenMax, lenGap, seqLong, TRANS, EMIS)
%
% Busca els trams de la sequencia d'entrada 'seqLong' amb màxima
% versemblaça, donat un model (matrius 'TRANS', 'EMIS') i les longituds
% minima 'lenMin' i maxima 'lenMax', així com la precissió de cerca
% 'lenGap'.
%
% Retorna una squencia 'seqLikelihoods' de la mateixa longitud que la
% sequencia d'entrada, amb la suma de les versemblances a cada posicio.
%

seqLikelihoods=zeros(size(seqLong));

len=length(seqLong);
for ind0=1:len-lenMin
    for l=lenMin:lenGap:lenMax
        ind1=ind0+l;
        if ind1<=len
            p = evaluateHMM(seqLong(ind0:ind1),TRANS, EMIS);
            for i=ind0:ind1;
                seqLikelihoods(i) = seqLikelihoods(i) + p;
            end
        end
    end
end

end

