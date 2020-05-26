function [new_interleaved, new_ts]=ca_TTL_interleaved(correct_folders)
%we want indexes of when chunks start and finish 
events=load('Event.mat');
events=events.Event;

%lets format these ttls. 
%first interleave since these are single pulses
interleaved = zeros(1,length(events)*2);
interleaved(1:2:end-1) = events(:,1);
interleaved(2:2:end)   = events(:,2);
interleaved=interleaved';
intervals=diff(interleaved);
common_interval=median(intervals(:,1));

%find recording chunks
jitter=0.0002;
raw_chunks=find(intervals>common_interval+jitter);

%remove single isolated TTLS
if raw_chunks(1)==1
   dr=diff(raw_chunks)==1;
   d=1;
   while dr(d)==1
         d=d+1;
   end
    starts=(d+1);
    raw_chunks(1:d)=[];

else
   starts=1; 
end

remove=[];
keep=[];
for i=1:length(raw_chunks)
    if (raw_chunks(i)+2)<=length(interleaved)
   
    if interleaved(raw_chunks(i)+1)>interleaved(raw_chunks(i))+common_interval+jitter ... 
        && interleaved(raw_chunks(i)+1)>interleaved(raw_chunks(i)+2)-common_interval-jitter
        keep=[keep; raw_chunks(i)+1]; %keep the stops
    end
    
    if interleaved(raw_chunks(i)+1)>interleaved(raw_chunks(i))+common_interval+jitter ...
            && interleaved(raw_chunks(i)-1)<interleaved(raw_chunks(i))-common_interval-jitter        
        remove=[remove; raw_chunks(i)];
        
    end
    end
end
stops=setdiff(raw_chunks,remove);
remove(diff(remove)==1)=[];
if length(remove)<length(stops)
    starts=[starts; keep];
else
    starts=[starts; remove+1];
end

%add last stop if needed
if starts(end)> raw_chunks(end)
    stops(end+1)=length(interleaved);
end
for s=1:length(starts)
    new_interleaved{s}=interleaved(starts(s):stops(s));
end
%%
%use timestamp.dat times
%read H folders only
c_dir=pwd;
cd ..
counter=1;
new_ts=[];
for f = 1:length(correct_folders)
    
    currD = correct_folders{f}; 
    cd(currD) 

    %for f=1 each folder
    timestamp=import_timestamp;
    ts=timestamp.sysClock;

    %ignore sysClocks not increasing
    go=1;        
    while go
        if ts(1)>ts(2)
            ts(1)=0;                
        else
            go=0;
        end
    end

    new_ts{counter}=ts;
    counter=counter+1;
    cd ..
    
end
cd (c_dir)
end








