function [outputArg1,outputArg2,depth_interp] = interp_vertical_cell(cellarray,depth_original,new_depth_interval)

% depth_original == vector with depth levels of original data
% depth_interp_interval== new interval of interpolation (ex . every 200 m)
% Extracts time series
% Mount back the original data into separate matrices
% Then interpolate vertically into new depth levels 
% Then builds a new cell array with the interpolated data, each cell containing a time series 


depth_interp=min(depth_original):new_depth_interval:max(depth_original); 
depth_interp=depth_interp';

PI_detrend_yearly=zeros(size(cellarray,1),size(cellarray,2),size(cellarray,3));
PI_detrend_interp=zeros(size(cellarray,1),size(cellarray,2),length(depth_interp));




for year=1:length(cellarray{150,150,1}) %this gives the length of the time series for a point I'm sure there is data
    year  
for l=1:size(cellarray,1)
    for c=1:size(cellarray,2)
       for  k=1:length(depth_original)
       
           if length((cellarray{l,c,k}))==30 %if there is a time series and not a NaN value
           PI_detrend_yearly(l,c,k)=cellarray{l,c,k}(year); % extracting the first point of the time series (year=1) for entire domain
           
           else
           PI_detrend_yearly(l,c,k)=NaN;
           end
       end
    
       dummy=squeeze(PI_detrend_yearly(l,c,:));
       dummy_interp=interpn(depth_original,squeeze(PI_detrend_yearly(l,c,:)),depth_interp);
       PI_detrend_interp(l,c,:)=dummy_interp;
    
       end
       
end

    
    strucdummy(year).data=PI_detrend_interp;
    strucdummy(year).Year=year;
    end
 



outputArg1 = strucdummy;


%%mount cell new series with interpolated data

for  l=1:size(cellarray,1)
     l
    for c=1:size(cellarray,2)
        for k=1:length(depth_interp)

            for year=1:length(cellarray{150,150,1})

                dummy2(year)=strucdummy(year).data(l,c,k);
            end

            PI_detrend_interp_series{l,c,k}=dummy2;

        end
    end
end



outputArg2 = PI_detrend_interp_series;
end

