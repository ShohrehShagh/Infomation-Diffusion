function [t]=Just_t(N,batch,term,L0,BL_t,prev_t)
t=zeros(N,1);
BL=max(size(term,2)-1,BL_t);
index=zeros(N,1);
if(batch==1)
    working_indices=L0+1:N;
else
    working_indices=find(isnan(prev_t));
end
for ii=1:length(working_indices)
    i=working_indices(ii);
    index(i)=find(rand <= cumsum(term(i,:))/sum(term(i,:)), 1);
end
t(index==BL+1)=NaN;
t(index~=(BL+1))=(batch-1)*BL_t+index(index~=(BL+1));
if(batch>1)
    t(~isnan(prev_t))=prev_t(~isnan(prev_t));
end
