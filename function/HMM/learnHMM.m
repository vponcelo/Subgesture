function [TRANS_EST, EMIS_EST] = learnHMM(nStates,seqs,hmmIters)
%
% Donada una sequencia 'seq' i el numero d'estats 'nStates', genera les
% etiquetes d'estat i apren les matrius de transferencia i emissio
% 
% 'seqs' no pot contenir cap 0 !

% Creem una matriu amb una fila per seqüència d'entrada, amb mida la mida
% de la seqüència més llarga
len=0;
if iscell(seqs)
    for i=1:length(seqs),
        if length(seqs{i})>len,
            len=length(seqs{i});
        end
    end
else
    len = length(seqs);
end    

seq=[];
uniformSeqs=[];
states=[];
stateCell={};

if iscell(seqs),
    for i=1:length(seqs),
        batch=seqs{i};
        [s1, s2]=size(batch);
        if(s2==1)
            batch=batch';
        end
        %seq=[seq batch];
        seq=[seq batch(1) batch(1) batch];

        %states=[states stateLabels(nStates,length(batch))];
        states=[states 1 1 stateLabels(nStates,length(batch))];

        uniformSeqs=[uniformSeqs;imresize(batch,[1 len],'nearest')];
        stateCell=[stateCell stateLabels(nStates,length(batch))];
    end
else
    batch=seqs;
    [s1, s2]=size(batch);
    if(s2==1)
        batch=batch';
    end
    %seq=[seq batch];
    seq=[seq batch(1) batch(1) batch];

    %states=[states stateLabels(nStates,length(batch))];
    states=[states 1 1 stateLabels(nStates,length(batch))];

    uniformSeqs=[uniformSeqs;imresize(batch,[1 len],'nearest')];
    stateCell=[stateCell stateLabels(nStates,length(batch))];
end

% seq ha de ser un vector fila:
%[s1, s2]=size(seq);
%if(s2==1)
%    seq=seq';
%end
%len=max(s1,s2);

% Es generen les etiquetes pels estats:
%states=stateLabels(nStates,len);
% S'estimen les matrius de transferencia i emissio
[TRANS_EST, EMIS_EST] = hmmestimate(seq, states);

% Es refina el resultat.
%[TRANS_EST, EMIS_EST] = hmmtrain(uniformSeqs, TRANS_EST, EMIS_EST,'Maxiterations',300','Verbose',true,'Tolerance',1e-6);
[TRANS_EST, EMIS_EST] = hmmtrain(seq, TRANS_EST, EMIS_EST,'Maxiterations',hmmIters','Verbose',false,'Tolerance',1e-6);

end

