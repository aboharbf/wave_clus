function parsed_chs = parse_data_NSx(filename,max_memo_GB,output_name,channels)
% Code only valid for recordings without pauses or data loss (Nsx.Data field can't be a cell)
% This code requires the file openNSx.m, from the NPMK, in the path. You can download the download NPMK from https://github.com/BlackrockMicrosystems/NPMK/releases
% max_memo_GB is an idea of the number of GB allocated for the data to be
% stored in RAM, so it is used to compute the number of segments in which
% the data should be split for processing

if ~exist('filename')
  filename = uigetfile('*.ns5');
end

if length(filename)<3 || (~strcmpi(filename(2:3),':\') && ...
    ~strcmpi(filename(1),'/') && ...
    ~strcmpi(filename(2),'/') && ...
    ~strcmpi(filename(1:2), '\\')&& ~strcmpi(filename(2:3),':/'))
  
  filename = [pwd filesep filename];
end

[filepath, file, fext] = fileparts(filename);
parsedFolder = [filepath '/' file '_parsed'];
fext = lower(fext(2:end));
nsx_ext = fext(end);
chext = ['.NC' nsx_ext];

NSx = openNSx(filename, 'report','noread');
nchan = NSx.MetaTags.ChannelCount;   % number of channels
sr = NSx.MetaTags.SamplingFreq;   % sampling rate
lts = sum(NSx.MetaTags.DataPoints);   % total data points
outfile_handles = cell(1,nchan); %some will be empty

if ~exist('output_name','var') || isempty(output_name)
  output_name = [parsedFolder '/' file '_NSX'];
end

% if ~isempty(outDirPath) MAYBE REDUNDENT NOW?
%   disp('Data already processed - returning appropriate paths')
%   outfile_path = fullfile({outDirPath.folder}, {outDirPath.name})';
%   return
% outDirPath = dir([parsedFolder '/'  file '_NSX_Ch*.NC5']);

% else

parsed_chs = [];
for i = 1:nchan
  c = NSx.MetaTags.ChannelID(i);
  if ismember(c,channels)
    parsed_chs(end+1) = c;
    outfile_handles{i} = fopen([output_name '_' num2str(c) chext],'w');
  end
end
DataPoints = NSx.MetaTags.DataPoints;

DataPoints = [0 DataPoints]; %added for use as starting index
samples_per_channel = ceil(max_memo/(nchan*length(NSx.MetaTags.DataPoints))/2);
for part = 2:length(DataPoints)
  N = DataPoints(part);   % total data points
  num_segments = ceil(N/samples_per_channel);
  fprintf('TimeStamps generated. Data will be processed in %d segments of %d samples each.\n',num_segments,min(samples_per_channel,N))
  for j=1:num_segments
    ini = (j-1)*samples_per_channel+1+DataPoints(part-1);
    fin = min(j*samples_per_channel,N)+DataPoints(part-1);
    tcum = tcum + toc;  % this is because openNSx has a tic at the beginning
    NSx = openNSx('read',filename,['t:' num2str(ini) ':' num2str(fin)]);
    for i = 1:nchan
      if ~isempty(outfile_handles{i})
        %BE AWARE THAT THE RAW DATA IN THE NS5 AND NC5 IS SCALED UP BY A FACTOR OF 4
        fwrite(outfile_handles{i},NSx.Data(i,:),'int16');
      end
    end
  end
  fprintf('Segment %d out of %d processed. Data Point Read = %d \n',j,num_segments,size(NSx.Data,2));
end
end
fclose('all');
tcum = tcum + toc;
fprintf('Total time spent in parsing the data was %s secs.\n',num2str(tcum, '%0.1f'));
end

fclose('all');
tcum = tcum + toc;
fprintf('Total time spent in parsing the data was %s secs.\n',num2str(tcum, '%0.1f'));

metadata_file = fullfile(pwd, [output_name '_TimeStamps.mat']);
metadata = struct;
if exist(metadata_file,'file')
  metadata = load(metadata_file);
end

if strcmp(nsx_ext,'5')
  if isfield(metadata,'parsed_chs')
    parsed_chs = union(parsed_chs,metadata.parsed_chs);
  end
  metadata.parsed_chs = parsed_chs;
  metadata.lts = lts;
  metadata.nchan = nchan;
  metadata.sr = sr;
  metadata.chext = chext;
  save(metadata_file, '-struct', 'metadata')
else
  if isfield(metadata,fext) && isfield(metadata.(fext),'sr')
    parsed_chs = union(parsed_chs,metadata.(fext).parsed_chs);
  end
  metadata.(fext).parsed_chs = parsed_chs;
  metadata.(fext).lts = lts;
  metadata.(fext).nchan = nchan;
  metadata.(fext).sr = sr;
  metadata.(fext).chext = chext;
  save(metadata_file, '-struct', 'metadata')
end
