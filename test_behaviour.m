%starting in animal direcctory
fs=dir();

for i=9:length(fs)    
    cd (fs(i).name)
    tic
    %1st classification of folders. Will have to do hypno for more careful analysis
    %2 is sw, 3 rem, 1 beh, 4nouse/further process
    %first take only folders with H
    ca_names = dir('H*');
    folders = {ca_names([ca_names.isdir]).name};
    c_folders=folders;
    %we have to order the folders
    for f=1: length(folders)
        times{1,f} = sscanf(folders{f},'H%d_M%d_S%d*');
    end
    wrong_digit=cellfun(@(x) x<10, times, 'UniformOutput', false);
    wrong_digit_folder=times(logical(cellfun(@sum, wrong_digit)));
    replacing_ind=find(cellfun(@sum, wrong_digit));

    for i=1: length (replacing_ind)
        if sum(wrong_digit{replacing_ind(i)})==1;
            if wrong_digit{replacing_ind(i)}(1)
                folders{replacing_ind(i)} = sprintf('H0%d_M%d_S%d',wrong_digit_folder{i}(1),wrong_digit_folder{i}(2),wrong_digit_folder{i}(3));
            elseif wrong_digit{replacing_ind(i)}(2)
                folders{replacing_ind(i)} = sprintf('H%d_M0%d_S%d',wrong_digit_folder{i}(1),wrong_digit_folder{i}(2),wrong_digit_folder{i}(3));
            else
                folders{replacing_ind(i)} = sprintf('H%d_M%d_S0%d',wrong_digit_folder{i}(1),wrong_digit_folder{i}(2),wrong_digit_folder{i}(3));
            end
        elseif sum(wrong_digit{replacing_ind(i)})==2;
             if wrong_digit{replacing_ind(i)}(1) && wrong_digit{replacing_ind(i)}(2)
                folders{replacing_ind(i)} = sprintf('H0%d_M0%d_S%d',wrong_digit_folder{i}(1),wrong_digit_folder{i}(2),wrong_digit_folder{i}(3));
             elseif wrong_digit{replacing_ind(i)}(1) && wrong_digit{replacing_ind(i)}(3)
                folders{replacing_ind(i)} = sprintf('H0%d_M%d_S0%d',wrong_digit_folder{i}(1),wrong_digit_folder{i}(2),wrong_digit_folder{i}(3));
             else
                folders{replacing_ind(i)} = sprintf('H%d_M0%d_S0%d',wrong_digit_folder{i}(1),wrong_digit_folder{i}(2),wrong_digit_folder{i}(3));
             end
        else
            folders{replacing_ind(i)} = sprintf('H0%d_M0%d_S0%d',wrong_digit_folder{i}(1),wrong_digit_folder{i}(2),wrong_digit_folder{i}(3));
        end
    end

    [~, correct_index]=sort(folders);
    correct_folders=c_folders(correct_index);
    folder_class=4+zeros(length(folders),1);
    for f=1: length(folders)
        if contains(correct_folders{f},'mix','IgnoreCase',true);
            folder_class(f)= 4;
        elseif contains(correct_folders{f},'swrem','IgnoreCase',true);
            folder_class(f)= 4;
        elseif contains(correct_folders{f},'remsw','IgnoreCase',true);
        folder_class(f)= 4;
        elseif contains(correct_folders{f},'rem','IgnoreCase',true);
           folder_class(f)= 3;
       elseif contains(correct_folders{f}, 'sw', 'IgnoreCase',true);
           folder_class(f)= 2;
       elseif contains(correct_folders{f}, 'l', 'IgnoreCase',true);
           folder_class(f)= 1;
       end
    end
    %%
    load('ms.mat');
    %let's look at the behav segment
    inx=find(folder_class==1);
    load('behav.mat');
    %ms to pass
    msB=[];
    %correct TS
    tf=sum(ms.timestamps0(1:inx-1));  
    msB.time=ms.time(tf+1:tf+ms.timestamps1(inx))-ms.time(tf+1); %timestamps0 for non-beh

    %slice frames
    msB.RawTraces=ms.RawTraces(tf+1:tf+ms.timestamps1(inx), :);%timestamps0 for non-beh
    msB.numFrames=tf+1+ms.timestamps1(inx)-(tf+1);%timestamps0 for non-beh
    msB.numNeurons=ms.numNeurons;
    msB= msExtractBinary(msB);
    %%
    %per cell loop
    try 
        cd ("spatial_analysis")
    catch
        mkdir spatial_analysis    
        cd ("spatial_analysis")
    end

    for c=1:100
%         LTA_bin{c} = msLinearTrackAnalysis_EVV(msB, behav, c,1,1);
%         name=sprintf('LTA_cell%d.fig', c);
%         saveas(gcf, name); close
%         LTA_raw{c} = msLinearTrackAnalysis_EVV(msB, behav, c,1, 2);
%         name=sprintf('LTA_raw_cell%d.fig', c);
%         saveas(gcf, name); close
        SF_bin{c} = msSpatialFiring_EVV_clean(msB, behav, c, 1,1);
%         name=sprintf('SF_cell%d.fig', c);
%         saveas(gcf, name); close
%         SF_raw{c} = msSpatialFiring_EVV_clean(msB, behav, c, 1,2);
%         name=sprintf('SF_raw_cell%d.fig', c);
%         saveas(gcf, name); close
    end
    
    for c=101:msB.numNeurons
         LTA_bin{c} = msLinearTrackAnalysis_EVV(msB, behav, c,0,1);
%         name=sprintf('LTA_cell%d.fig', c);
%         saveas(gcf, name); close
        LTA_raw{c} = msLinearTrackAnalysis_EVV(msB, behav, c,0, 2);
%         name=sprintf('LTA_raw_cell%d.fig', c);
%         saveas(gcf, name); close
        SF_bin{c} = msSpatialFiring_EVV_clean(msB, behav, c, 0,1);
%         name=sprintf('SF_cell%d.fig', c);
%         saveas(gcf, name); close
        SF_raw{c} = msSpatialFiring_EVV_clean(msB, behav, c, 0,2);
%         name=sprintf('SF_raw_cell%d.fig', c);
%         saveas(gcf, name); close
    end
    save("spatial_analysis.mat", 'LTA_bin', 'LTA_raw', 'SF_bin', 'SF_raw');
    cd ..
    toc
    cd ..
    clearvars -except i fs
end
%%
isPlaceCell= cellfun(@(c) c.IsPlaceCell, SF_bin, 'UniformOutput', false);
isPlaceCell=cell2mat(isPlaceCell);
sum(isPlaceCell)/length(isPlaceCell)

PlaceFieldStability= cellfun(@(c) c.PlaceFieldStability, SF_bin, 'UniformOutput', false);
PlaceFieldStability= cellfun(@mean,  PlaceFieldStability, 'UniformOutput', false);
PlaceFieldStability_PC=cell2mat(PlaceFieldStability(isPlaceCell));
mean(PlaceFieldStability_PC)

CellFiringProb  = cellfun(@(c) c.CellFiringProb, SF_bin, 'UniformOutput', false);
CellFiringProb=cell2mat(CellFiringProb(isPlaceCell));
mean(CellFiringProb)

numPlaceFields = cellfun(@(c) c.numPlaceFields, SF_bin, 'UniformOutput', false);
numPlaceFields=cell2mat(numPlaceFields);


PlaceFieldArea  = cellfun(@(c) c.PlaceFieldArea, SF_bin, 'UniformOutput', false);
PlaceFieldArea=(PlaceFieldArea(numPlaceFields==1));
PlaceFieldArea=cell2mat(vertcat(PlaceFieldArea{:}));
mean(PlaceFieldArea);

InFieldActivityRatio  = cellfun(@(c) c.InFieldActivityRatio, SF_bin, 'UniformOutput', false);
InFieldActivityRatio=cell2mat(InFieldActivityRatio(numPlaceFields==1));
mean(InFieldActivityRatio);
