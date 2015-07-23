function W = genRepresentation(dataw,CVAL)
    Dataw =  zeros(size(dataw,1)-1,size(dataw,2));
    for j = 2 : size(dataw,1)                
        Dataw(j-1,:) = mean(dataw(1:j,:));
    end            
    temp=Dataw';
% %     %temp = vl_homkermap(Dataw',2,'kchi2');                
%     temp = (sqrt(abs(Dataw')));
    W_fow = liblinearsvr(temp',CVAL,2); clear temp;				
    order = 1:size(dataw,1);
    [~,order] = sort(order,'descend');
    dataw = dataw(order,:);
    Dataw =  zeros(size(dataw,1)-1,size(dataw,2));
    for j = 2 : size(dataw,1)                
        Dataw(j-1,:) = mean(dataw(1:j,:));
    end            
%     %temp = vl_homkermap(Dataw',2,'kchi2');                
%     temp = sqrt(abs(Dataw'));
    temp=Dataw';
    W_rev = liblinearsvr(temp',CVAL,2); clear temp;				              
    W = [W_fow ; W_rev];  

end
