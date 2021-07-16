
function [delta, year1_field,year2_field] = do_delta(varname,mode,varargin)

%in the future do this for each in the sequence ear so I can do an
%animation, thus loading the field only once into memory and reiterating
%over all 

% Do Delta between first and last year of simulation
% Optional arguments 'year1' and 'year2'
% delta = do_delta('templvl') % delta between last year of run and year1
% default value for year1 = 1
% default value for year2 = 481 (length of NorESM1-ME ramp simulation)
% delta = do_delta('templvl','year1',1, 'year2',140) %% delta between year 140 and year 1
   
field = load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Binded_var_files/NorESM_filtered_%s',varname));
var_field = field.var_field;
clear field

% unpacking to see if it works

%first get exemple of timeseries length and other dimensions
for line = 1:size(var_field,1)
    for col = 1:size(var_field,2)
          for k= 1:size(var_field,3)
              if ~isnan(var_field{line,col,k})
                  time_series_length=length(var_field{line,col,k});
                  break
              end
          end
    end
end


defaultyear1 = 1;
defaultyear2 = time_series_length;

p = inputParser;
addRequired(p,'varname')
addRequired(p,'mode')
addParameter(p,'year1', defaultyear1)
addParameter(p,'year2', defaultyear2)

parse(p,varname,mode, varargin{:});

%Extract values from parser
year1 = p.Results.year1 ;
year2 = p.Results.year2 ;
mode = p.Results.mode;

extract_y1=NaN(size(var_field,1),size(var_field,2),size(var_field,3));
extract_y2=NaN(size(var_field,1),size(var_field,2),size(var_field,3));


if strcmp(mode,'year_pairs')
    sprintf('Mode is %s', mode)

for line = 1:size(var_field,1)
    for col = 1:size(var_field,2)
          for k= 1:size(var_field,3)
              if length(var_field{line,col,k}) == time_series_length
                  extract_y1(line,col,k) = var_field{line,col,k}(year1);
                  extract_y2(line,col,k) = var_field{line,col,k}(year2);
              else
                extract_y1(line,col,k) = NaN;
                extract_y2(line,col,k) = NaN;

              end
          end
    end
end

year1_field = extract_y1;
year2_field = extract_y2;
delta = extract_y2 - extract_y1;

end


if strcmp(mode,'complete')
    sprintf('Mode is %s', mode)
   
    dummy = cell(size(var_field,1),size(var_field,2),size(var_field,3));
   
    for line = 1:size(var_field,1)
        for col = 1:size(var_field,2)
              for k = 1:size(var_field,3)
                  if length(var_field{line,col,k}) == time_series_length
                      
                      first_elemen = var_field{line,col,k}(year1);
                      extract_year1 = repmat(first_elemen,1,time_series_length);
                      
                      extract_y2 = var_field{line,col,k}(:)'; %all years
                      dummy{line,col,k} = extract_y2 - extract_year1; %vector series of deltas
                      
                  else
                    dummy{line,col,k} = NaN;

                  end
              end
        end
    end
    
    year1_field = [];
    year2_field = [];
    delta = dummy;

end




              
              
              

