function [ Regressor, Vol4, Vol3, Vol2, Vol1, PPI, AudHab ] = ASD_standard_regressor()
%%
%
%
%
%
%
%   Assumptions:
%       Length of imaging is 1510
%
%


stim_num_frames = 1510;
%Building the regressor by ID the stim times and then applying to data
spike=[0,1.69644104899772,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
spike=spike/max(spike);
Regressor=zeros(1, stim_num_frames);

%classify timepoints of stimulus - Vol 4 is the loudest Vol 1 is the
%quietist, then there is PPI and Habituation

Vol4=[200, 400, 520, 840, 1080]; %-6dB
for i=Vol4
    Regressor(1,i-1:i-1+length(spike)-1)=spike';
end
Vol3=[120, 440, 680, 920, 1040]; %-12dB
for i=Vol3
    Regressor(1,i-1:i-1+length(spike)-1)=spike';
end
Vol2=[160, 320, 600, 720, 960]; %-18dB
for i=Vol2
    Regressor(1,i-1:i-1+length(spike)-1)=spike';
end
Vol1=[280, 560, 760, 880, 1000]; %-24dB
for i=Vol1
    Regressor(1,i-1:i-1+length(spike)-1)=spike';
end
PPI=[240, 360, 480, 640, 800];
for i=PPI
    Regressor(2,i-1:i-1+length(spike)-1)=spike';
end

AudHab=(10 * 60 * 2 : 6 : 11 * 60 * 2 - 6); % Only includes first habituation block
second_block = (12 * 60 * 2 : 6 : 12.5 * 60 * 2 -6);
AudHab = [AudHab, second_block];
for i=1:numel(AudHab)
    Regressor(3,AudHab(i):AudHab(i)+length(spike)-1)=spike';
end

Regressor = Regressor(:, 1:1510);

end