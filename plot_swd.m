frame=20714;%+83
[ms] = msExtractBinary(ms);
[~, ind]=sort(ms.RawTraces(frame,:));
[~, ind2]=sort(ms.Binary(frame,:));

figure;pcolor((ms.RawTraces(frame-10:frame+19,ind))');shading interp;
caxis([-0.2 0.2])
figure;imagesc((ms.RawTraces(frame-10:frame+19,ind))');

mframe=mean(ms.RawTraces,2);
figure; plot(mframe(frame-10:frame+19));
