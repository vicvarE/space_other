
%initialize
[new_interleaved, new_ts]=ca_TTL_interleaved();
num_recs=length(new_ts);
load('swdH.mat');
load('ms.mat');
load('Downsampled_Sampling_Frequency.mat');
Fs=Target_Sampling_Frequency;
channel= sprintf('CSC00%d.mat', swdH.channel);
load(channel) ;
resampled_data_z=zscore(resampled_data_mV);

[ms] = msExtractBinary(ms);

ztraces=zscore(ms.RawTraces);
mean_cells=(mean(ztraces,2));
cam_Fs=30;
if num_recs ~= length(new_interleaved)
    warning("DAQ stopped sending signals, you need to merge certain TTLs")
    return
end

%%
%Analysis
accum_reclen=0;

swd_ind=[];
swd_peak=[]; 
swd_peakz=[]; 
pop_peak=[]; 
ct_p=[];
ct_nm=[];
for r=1:num_recs
    interleaved=new_interleaved{r};
    swd_rec=swdH.pksxs>interleaved(1) & swdH.pksxs<interleaved(end); 
    swds=swdH.pksxs(swd_rec); %skip recs with no swds
    if isempty(swds)
        continue
    end
    ts=new_ts{r};              
    swd_indices=[];
    swd_peakTEMP=[]; 
    swd_peakzTEMP=[];
    pop_peakTEMP=[];
    cells_timep=[];
    cells_timenm=[];
    
    n_dir =(sprintf('rec%d',r));
    if ~exist(n_dir, 'dir')
        mkdir(n_dir)
    end
    cd (n_dir)
    for s=1: length(swds)

        corrected_swd=(swds(s)-interleaved(1));
        index1=accum_reclen+round((swds(s)-interleaved(1))*cam_Fs); %have to sum the lenght of previous recording
        %[~, index2] = min(abs(ts/1000-corrected_swd)); %closest frame
        diff_frame=(ts/1000-corrected_swd);
        diff_frame(diff_frame < 0 ) = NaN;
        [~, index2] = min(diff_frame); %closest largest frame
        index2=index2+accum_reclen;
        %[~, index3] = min(abs(interleaved-swds(i))); %closest ttl
        difference=(interleaved-swds(s));
        difference(difference < 0 ) = NaN;
        [~, index3] = min(difference); %closest largest ttl     
        index3=index3+accum_reclen;
       
        figure;
        subplot(3,3,1);plot(resampled_data_mV(round(swds(s)*Fs)-Fs/2:round(swds(s)*Fs)+(Fs/2)-1));
        maxp=max(resampled_data_mV(round(swds(s)*Fs)-Fs/4:round(swds(s)*Fs)+(Fs/4)-1));
        maxpz=max(resampled_data_z(round(swds(s)*Fs)-Fs/4:round(swds(s)*Fs)+(Fs/4)-1));
        
        [maxpop,index4]=max(mean_cells(index2-cam_Fs/2:index2+(cam_Fs/2)-1));
        index4=index4+index2-((cam_Fs/2)+1);
        subplot(3,3,2);plot(mean_cells(index4-cam_Fs/2:index4+(cam_Fs/2)-1));

        [~, ind]=sort(ms.Binary(index4,:));
        subplot(3,3,4);pcolor((ms.Binary(index4-cam_Fs/2:index4+(cam_Fs/2)-1,ind))');shading interp;
        subplot(3,3,7); PlotCellsColorIndexv2(ms.SFPs, ind, ms.Binary(index4,:)); 

        pos=sum(ms.Binary(index4,:)==1);
        nm=sum(ms.Binary(index4,:)==0);
        subplot(3,3,3); pie([pos nm]);legend({'pos', 'nm'});
        
        posb=sum(sum(ms.Binary(index4-5:index4-1,:)==1)>0);
        nmb=sum(sum(ms.Binary(index4-5:index4-1,:)==0)>0);
        subplot(3,3,6); pie([posb nmb]);legend({'posb', 'nmb'})
        
        cells_time_pos=sum(ms.Binary(index4-cam_Fs/2:index4+(cam_Fs/2)-1,:)==1,2);
        subplot(3,3,5); bar(cells_time_pos);        
        
        cells_time_nm=sum(ms.Binary(index4-cam_Fs/2:index4+(cam_Fs/2)-1,:)==0,2);
        subplot(3,3,8); bar(cells_time_nm);
        
        drawnow;
        swd_indices(s, :)=[index1 index2 index3 index4];        
        swd_peakTEMP(s,:)=maxp;
        swd_peakzTEMP(s,:)=maxpz;
        pop_peakTEMP(s,:)=maxpop;
        cells_timep(:,s)=cells_time_pos;
        cells_timenm(:,s)=cells_time_nm;        
        
         n_fig =(sprintf('swd%d.fig',s));
         savefig(n_fig);
    end
    cd ..
    accum_reclen=accum_reclen+length(interleaved);
    swd_ind{r}=swd_indices;   
    swd_peak{r}=swd_peakTEMP; 
    swd_peakz{r}=swd_peakzTEMP; 
    pop_peak{r}=pop_peakTEMP;
    ct_p{r}=cells_timep;
    ct_nm{r}=cells_timenm;
    close all
end
%%
%if we dont care about state merge cells
ct_pm=cat(2, ct_p{:});
figure; bar((mean(ct_pm,2)/ms.numNeurons)*100);

peakz=cat(1,swd_peakz{:});
coef = polyfit(repmat(ct_pm(15,:), size(peakz',1), 1), peakz', 1);
figure;scatter(ct_pm(15,:), peakz); 
h = refline(coef(1), coef(2));
[r,p] = corrcoef(ct_pm(15,:),peakz');

p_peak=cat(1,pop_peak{:});
coef = polyfit(repmat(ct_pm(15,:), size(p_peak',1), 1), p_peak', 1);
figure;scatter(ct_pm(15,:), p_peak); 
h = refline(coef(1), coef(2));
[r,p] = corrcoef(ct_pm(15,:),p_peak');

coef = polyfit(repmat(peakz', size(p_peak',1), 1), p_peak', 1);
figure;scatter(peakz, p_peak); 
h = refline(coef(1), coef(2));
[r,p] = corrcoef(peakz,p_peak);


index=cat(1, swd_ind{:});
index=index(:,4);
for i=1:length(index)
    bin_swd_matrix(:,:,i)=ms.Binary(index(i)-cam_Fs/2:index(i)+(cam_Fs/2)-1,:);    
end
a=bin_swd_matrix(:,1,:);
a=permute(a,[1,3,2]);
figure; imagesc(a');

p_firing=mean(bin_swd_matrix,3);
figure; bar(p_firing(:,1));

[s_pf, ins]=sortrows(p_firing',15, 'descend');
figure; imagesc(s_pf);
caxis([0,1])
ylim([0 400])

b=a(14,:);
c=1:1:113;
d=c(b==1);
figure;histogram(d);

%%
%plot mean trace of all neurons + background