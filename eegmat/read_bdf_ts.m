function [ts, channel_names] = read_bdf_ts(bdffile, channels, start, stop)

if exist(bdffile) ~= 2
    error('BDF file not found');
end

hdr = read_bdf_hdr(bdffile);

%% Parse Inputs
switch nargin
    case 1
        channels = 'all';
        start = 0;
        stop = hdr.n_records - 1;
    case 2
        start = 0;
        stop = hdr.n_records - 1;
    case 3
        stop = hdr.n_records - 1;
    otherwise
end

if ischar(channels)
    switch channels
        case 'all'
            channels = 1:hdr.n_channels;
        case 'eeg'
            channels = find(strcmp(hdr.reserved,'MON'));
        case 'misc'
            channels = setdiff(1:hdr.n_channels, ...
                [find(strcmp(hdr.reserved,'MON')), ...
                find(strcmp(hdr.reserved,'TRI'))]);
        case 'trig'
            channels = find(strcmp(hdr.reserved,'TRI'));
        otherwise
            channels = find(strcmp(hdr.reserved,channels));
    end
    if isempty(channels)
        error('Cannot parse channels');
    end
else
    if length(setdiff(channels, 1:hdr.n_channels)) ~= 0
        error('Cannot parse channels.');
    end
end

if start > hdr.n_records - 1
    error('Invalid start time');
end
if stop > hdr.n_records || stop <= start
    error('Invalid stop time');
end

%% Initialize BDF file
fid = fopen(bdffile, 'r');

%% Read data
sample_rate = unique(hdr.n_samples(channels));
if numel(sample_rate) > 1
    error('Some of the specified channels have different sample rates, please read these in separately');
end

recs_read = floor(start):floor(stop);
n_recs = numel(recs_read);

ts = zeros(sample_rate*n_recs, numel(channels));
channel_names = hdr.channel_labels(channels);

for i = 1:n_recs
    for j = 1:numel(channels)
        offset = hdr.hdr_bytes + sum(hdr.n_samples)*recs_read(i)*3 + ...
            (channels(j)-1)*sample_rate*3;
        fseek(fid, offset, 'bof');
        temp = fread(fid, 3*sample_rate);
        temp = reshape(temp, 3, [])';
        if channels(j) == hdr.n_channels; temp(:,end) = 0; end
        temp = sum(bsxfun(@times, temp, [1, 256, 256^2]),2);
        ineg = temp >= 256*256*128;
        temp(ineg) = temp(ineg) - 256^3;
        
        if channels(j) == hdr.n_channels
            gain = 1;
        else
            gain = ((hdr.phys_max(channels(j))-hdr.phys_min(channels(j))) ./ ...
                (hdr.digi_max(channels(j))-hdr.digi_min(channels(j))))/1000000;
        end
        temp = temp * gain;
        ts(((i-1)*sample_rate+1):(i*sample_rate),j) = temp;
    end
end

trim1 = round((start - floor(start)) * sample_rate)+1;
trim2 = sample_rate - round((stop-floor(stop)) * sample_rate) + 1;
ts = ts(trim1:(size(ts,1)-trim2), :);

fclose(fid);
