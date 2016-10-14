clc
clear all
for r=1%:10
    for k=3:4
        load(sprintf('Samples/Samples_online_tza_%d_4_%d',k,r))
        load(sprintf('input_datas/Data_%d_%d.mat',k,r))
        for i=1:size(Samples_To_Save.parents,2)
            zs=Samples_To_Save.parents(:,i);
            zs(isnan(zs))=inf;
            [a_z,~,c_z] = unique(zs,'rows');
            percentage1=accumarray(c_z,1)/size(Samples_To_Save.parents,1)*100;
            if(isempty(find(a_z(:,1)==Real.z(i))))
                percentage_correct(i)=0;
                k
            else
                percentage_correct(i)=percentage1(find(a_z(:,1)==Real.z(i)));
            end
        end
        percentage_correct
        mean_percentage_correct(r,k)=mean(percentage_correct);
        %
        load(sprintf('Samples/Samples_online_za_%d_4_%d',k,r))
        load(sprintf('input_datas/Data_%d_%d.mat',k,r))
        for i=1:size(Samples_To_Save.parents,2)
            zs=Samples_To_Save.parents(:,i);
            zs(isnan(zs))=inf;
            [a_z,~,c_z] = unique(zs,'rows');
            percentage1=accumarray(c_z,1)/size(Samples_To_Save.parents,1)*100;
            if(isempty(find(a_z(:,1)==Real.z(i))))
                percentage_correct(i)=0;
                k
            else
                percentage_correct(i)=percentage1(find(a_z(:,1)==Real.z(i)));
            end
        end
        percentage_correct
        mean_percentage_correct(r,k)=mean(percentage_correct);
        %
    end
end
%mean(mean_percentage_correct)