
events=load('Event.mat');
events=events.Event;

%lets format these ttls
%first interleave since these are single pulses
interleaved = zeros(1,length(events)*2);
interleaved(1:2:end-1) = events(:,1);
interleaved(2:2:end)   = events(:,2);
interleaved=interleaved';
intervals=diff(interleaved);
common_interval=median(intervals(:,1));

%find recording chunks
jitter=0.0001;
raw_chunks=find(intervals>common_interval+jitter);

%remove single isolated TTLS
if raw_chunks(1)==1
   raw_chunks(1)=[];
   starts=2;
else
   starts=1; 
end

remove=[];
for i=1:length(raw_chunks)
    if interleaved(raw_chunks(i)+1)>interleaved(raw_chunks(i))+common_interval+jitter ...
            && interleaved(raw_chunks(i)-1)<interleaved(raw_chunks(i))-common_interval-jitter        
        remove=[remove; raw_chunks(i)];
    end
end
stops=setdiff(raw_chunks,remove);
remove(diff(remove)==1)=[];
starts=[starts; remove+1];

%add last stop if needed
if starts(end)> raw_chunks(end)
    stops(end+1)=length(interleaved);
end
%%
%use timestamp.dat times
D = dir;
D=D([D.isdir]);
for f = 4:length(D)
    currD = D(f).name; 
    cd(currD) 
  
    %for f=1 each folder
    timestamp=import_timestamp;
    ts=timestamp.sysClock;
    interval_ts=diff(ts);
    %ignore sysClocks not increasing
    go=1;
    frameswo_clock=0;
    while go
        if ts(1)>ts(2)
            ts(1)=[];
            frameswo_clock=frameswo_clock+1;
        else
            go=0;
        end
    end
    ff=f-3;
    new_ts=interleaved(starts(ff)); %in s
    new_ts=[new_ts; new_ts(1)+(ts/1000)]; %orig in ms
    missing_frames=round((interleaved(stops(ff))-new_ts(end))/0.003);
    if ts(end)/1000<(interleaved(stops(ff))-interleaved(starts(ff)))
        msg = ["last ", num2str(missing_frames), "frame(s) is/are missing"];
        disp(msg)
    end
    test=[new_ts,interleaved(starts(ff):stops(ff))];
    test2=diff(test);
    comm_delay=mean(test(:,1)-test(:,2));
    cd ..
end
        
%%        
%calculate lengths of recordings to check with ms
rec_len=(stops)-(starts)+1;
if sum(rec_len==ms.timestamps)~=length(rec_len)
    %theres an offset in some/all recordings
    whichrecs=find(rec_len~=ms.timestamps);
    howmuch= rec_len(whichrecs)-ms.timestamps(whichrecs);
    
    for r=1:length(whichrecs)
    %import table
        checkframes=all(diff(timestamp.frameNum(2:end))==1);
        checkbuffer=all(diff(timestamp.buffer(2:end))==0);
        if checkframes==1 && checkbuffer==1        
            frame_interval=diff(timestamp.sysClock(3:end));
            checkclock=(diff(frame_interval)<30) & (diff(frame_interval)>-30);%skipped frame
            skip_ind=find(checkclock~=1);
            for s=1:length(skip_ind)
                messed(:,s)=frame_interval(skip_ind(s)-1:skip_ind(s)+1);
            end
        else 
            warning("problem with recording buffer");
            break
        end
    end
end

%find value pairs
stops=raw_chunks(diff(raw_chunks)==1);
starts=stops+1;
if raw_chunks(end)== starts(end)
    stops(end+1)=length(events)*2;    
else
    warning("deal with end of ttls, there may be extra solitaires")
end
if raw_chunks(1)==1 && raw_chunks(2)==stops(1)
    starts=[length(events)+1; starts];
else
    warning("deal with beginning of ttls, there may be extra solitaires")
end








