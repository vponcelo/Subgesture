function W = genRepresentation(data,CVAL)
    if nargin < 2
	CVAL = 1;
    end	
    OneToN = [1:size(data,1)]';    
    Data = cumsum(data);
    Data = Data ./ repmat(OneToN,1,size(Data,2));
    W_fow = liblinearsvr(getNonLinearity(Data),CVAL,2); clear Data; 			
    order = 1:size(data,1);
    [~,order] = sort(order,'descend');
    data = data(order,:);
    Data = cumsum(data);
    Data = Data ./ repmat(OneToN,1,size(Data,2));
    W_rev = liblinearsvr(getNonLinearity(Data),CVAL,2); 			              
    W = [W_fow ; W_rev]; 
end

function Data = getNonLinearity(Data)
%     Data = sign(Data).*sqrt(abs(Data));
%     Data = sign(Data).*sqrt(Data);
%     Data = sign(Data).*Data;
    %Data = vl_homkermap(Data',2,'kchi2');    
end

% function W = genRepresentation(dataw,CVAL)
%     Dataw =  zeros(size(dataw,1)-1,size(dataw,2));
%     for j = 2 : size(dataw,1)                
%         Dataw(j-1,:) = mean(dataw(1:j,:));
%     end            
%     temp=Dataw';
% % %     %temp = vl_homkermap(Dataw',2,'kchi2');                
% % %     temp = (sqrt(abs(Dataw')));
% % %     temp = sqrt(Dataw');
%     W_fow = liblinearsvr(temp',CVAL,2); clear temp;				
%     order = 1:size(dataw,1);
%     [~,order] = sort(order,'descend');
%     dataw = dataw(order,:);
%     Dataw =  zeros(size(dataw,1)-1,size(dataw,2));
%     for j = 2 : size(dataw,1)                
%         Dataw(j-1,:) = mean(dataw(1:j,:));
%     end            
% %     %temp = vl_homkermap(Dataw',2,'kchi2');
%     temp=Dataw';
% % %     temp = sqrt(abs(Dataw'));
% % %     temp = sqrt(Dataw');
%     W_rev = liblinearsvr(temp',CVAL,2); clear temp;				              
%     W = [W_fow ; W_rev];
% end
