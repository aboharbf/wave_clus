function [outfile_path] = parse_data_NSx(filename,max_memo_GB)
% Code only valid for recordings without pauses or data loss (Nsx.Data field can't be a cell)
% This code requires the file openNSx.m, from the NPMK, in the path. You can download the download NPMK from https://github.com/BlackrockMicrosystems/NPMK/releases
% max_memo_GB is an idea of the number of GB allocated for the data to be
% stored in RAM, so it is used to compute the number of segments in which
% the data should be split for processing

%% Process filename
if ~exist('filename','var')
  [fileN, path] = uigetfile('*.ns5');
  filename = [path filesep fileN];
end

if length(filename)<3 || (~strcmpi(filename(2:3),':\') && ...
    ~strcmpi(filename(2:3),':/') && ...
    ~strcmpi(filename(1),'/') && ...
    ~strcmpi(filename(1:2), '\\'))
  
  filename= [pwd filesep filename];
end


[filepath, file, ~] = fileparts(filename);
parsedFolder = [filepath '/' file '_parsed'];
outDirPath = dir([parsedFolder '/'  file '_NSX_Ch*.NC5']);

%% If the file is already parsed, return this path, otherwise, parse the file.
if ~isempty(outDirPath)
  disp('Data already processed - returning appropriate paths')
  for i = 1:length(outDirPath)
    outfile_path{i} = [outDirPath(i).folder filesep outDirPath(i).name];
  end
  return
else
  with_memory=true;
  try
    memory;
  catch
    with_memory=false;
  end
  if with_memory
    [userview,systemview] = memory;
    memo_avaible = floor(systemview.PhysicalMemory.Available*0.80);
    if exist('max_memo_GB','var') && ~isempty(max_memo_GB)
      max_memo = max_memo_GB*(1024)^3;
      if max_memo > memo_avaible
        error('max_memo_GB > 80% of Physical Memory Available')
      end
    else
      max_memo = memo_avaible;
    end
  else
    max_memo = max_memo_GB*(1024)^3;
  end
  
  tcum=0;
  
  mkdir(parsedFolder);
  NSx = openNSx(filename, 'report','noread');
  nchan = NSx.MetaTags.ChannelCount;   % number of channels
  sr = NSx.MetaTags.SamplingFreq;   % sampling rate
  lts = NSx.MetaTags.DataPoints;   % total data points
  % NSx.MetaTags.ChannelID;   % ChannelID
  % NSx.MetaTags.FileSpec   % Version
  % NSx.MetaTags.DataDurationSec   % total length in secs
  
  samples_per_channel = ceil(max_memo/nchan/2);
  num_segments = ceil(lts/samples_per_channel);
  
  for i = 1:nchan
    outfile_path{i} = [parsedFolder '/'  file '_NSX_Ch' num2str(NSx.MetaTags.ChannelID(i)) '.NC5'];
    outfile_handles{i} = fopen(outfile_path{i},'w');
  end
  
  TimeStamps=linspace(0,(lts-1)*1e6/sr,lts); %TimeStamps in microsec, with 0 corresponding to the first sample
  save([parsedFolder filesep file '_NSX_TimeStamps'],'TimeStamps','lts','nchan','sr');
  clear TimeStamps;
  fprintf('TimeStamps generated. Data will be processed in %d segments of %d samples each.\n',num_segments,min(samples_per_channel,lts))
  
  for j=1:num_segments
    ini = (j-1)*samples_per_channel+1;
    fin = min(j*samples_per_channel,lts);
    tcum = tcum + toc;  % this is because openNSx has a tic at the beginning
    NSx = openNSx('read',filename,['t:' num2str(ini) ':' num2str(fin)]);
    for i = 1:nchan
      fwrite(outfile_handles{i},NSx.Data(i,:),'int16');
    end
    fprintf('Segment %d out of %d processed. Data Point Read = %d \n',j,num_segments,size(NSx.Data,2));
  end
  
  fclose('all');
  tcum = tcum + toc;
  fprintf('Total time spent in parsing the data was %s secs.\n',num2str(tcum, '%0.1f'));
end


%fprintf('BE AWARE THAT THE RAW DATA IN THE NS5 AND NC5 IS SCALED UP BY A FACTOR OF 4.\n');

