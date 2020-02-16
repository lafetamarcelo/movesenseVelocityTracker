classdef sensor < handle
    %SENSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess=private)
        % Boolean flags for control
        resampled_flag               % Flag to control resampled data
        digital_data_flag            % Flag to control if data is digital
    end
    
    properties
        Name           % Name of the sensor used
        Timeseries     % The time series data collected
        SampleRate     % The sample ratge of the sensor (Hz)
        sampleTime     % The sample time of the sensor with variation
        dataSize       % Size of the data collected
        History        % The time series manipulation history
    end
    
    methods
        function obj = sensor(name,filename,multiplesensors)
            % This is the object constructor. It gets a file direction and
            % generates the sensor structure from the read data file
            % method.
            
            obj.Name = name; % Define the name
            if (multiplesensors) 
               pattern = '%s %s %s %s %s %s %s %s';
            else
               pattern = '%s %s %s %s %s %s %s';
            end
            % Read the csv file and get the sensor data
            fileid = fopen(filename);
            data = textscan(fileid, pattern,...
                'Delimiter',',', ...
                'EmptyValue',0);
            fclose(fileid);
            % Create the data structure
            obj.Timeseries.x = 0; obj.Timeseries.y = 0;
            obj.Timeseries.z = 0; obj.Timeseries.t = 0;
            % Get the info from raw data
            raw_data_size = length(data{2}) - 1;
            if (multiplesensors)
                % If there are multiple sensors
                obj.dataSize = 1;
                for i = 2 : (raw_data_size - 1)
                   if strcmp(data{1}(i),obj.Name)
                      obj.Timeseries.t(obj.dataSize) = ...
                          str2double(string(data{2}(i))); 
                      obj.Timeseries.x(obj.dataSize) = ...
                          str2double(string(data{3}(i)) + ...
                          "." + string(data{4}(i))); 
                      obj.Timeseries.y(obj.dataSize) = ...
                          str2double(string(data{5}(i)) + ...
                          "." + string(data{6}(i))); 
                      obj.Timeseries.z(obj.dataSize) = ...
                          str2double(string(data{7}(i)) + ...
                          "." + string(data{8}(i)));
                      obj.dataSize = obj.dataSize + 1;
                   end
                end
            else
                % If there is only one sensor
                for i = 2 : (raw_data_size - 1)
                  obj.Timeseries.t(obj.dataSize) = ...
                      str2double(string(data{1}(i))); 
                  obj.Timeseries.x(obj.dataSize) = ...
                      str2double(string(data{2}(i)) + ...
                      "." + string(data{3}(i))); 
                  obj.Timeseries.y(obj.dataSize) = ...
                      str2double(string(data{4}(i)) + ...
                      "." + string(data{5}(i))); 
                  obj.Timeseries.z(obj.dataSize) = ...
                      str2double(string(data{6}(i)) + ...
                      "." + string(data{7}(i)));
                  obj.dataSize = obj.dataSize + 1;
                end
            end
            
            % Compute the properties of the sensor
            obj.sampleTime.Value = mean(diff(obj.Timeseries.t));
            obj.sampleTime.Variation = std(diff(obj.Timeseries.t));
            obj.sampleTime.unity = "milliseconds";
            obj.digital_data_flag = (obj.sampleTime.Variation < 1E-06);
            
            % Show the sensor information logging
            obj.summary()
        end
        
        % Preprocessing of sensor
        function resample_series(obj, varargin)
           %
           %
           
           cfg.Enhance = 1;
           if ~isempty(varargin)
              for k = 1 : 2 : length(varargin)
                  cfg.(string(varargin{k})) = varargin{k+1};
              end
           end
           
           desired_freq = cfg.Enhance / obj.sampleTime.Value;
           
           time_data = obj.Timeseries.t;
           axis_data = [obj.Timeseries.x; ...
                        obj.Timeseries.y; ...
                        obj.Timeseries.z]';
           
           [data, time] = resample(axis_data, time_data, desired_freq);
           obj.digital_data_flag = true;
           obj.resampled_flag = true;
           
           % Save history of data
           obj.saveState("Resample");
           
           % Reatribute the time series value
           obj.Timeseries.t = time;
           obj.Timeseries.x = data(1:end,1)';
           obj.Timeseries.y = data(1:end,2)';
           obj.Timeseries.z = data(1:end,3)';
           obj.sampleTime.Value = mean(diff(time));
           obj.sampleTime.Variation = std(diff(time));
           obj.SampleRate = 1000 / obj.sampleTime.Value;
           obj.dataSize = length(obj.Timeseries.t);
           
           % Plot the new series summary
           obj.summary();
           
        end
        function denoise_series(obj, varargin)
            %
            %
            cfg = struct();
            if ~isempty(varargin)
              for k = 1 : 2 : length(varargin)
                  cfg.(string(varargin{k})) = varargin{k+1};
              end
            end
            
            % Save current state
            obj.saveState("Denoise");
            
            % Denoise each axis signal
            obj.Timeseries.x = wdenoise(obj.Timeseries.x);
            obj.Timeseries.y = wdenoise(obj.Timeseries.y);
            obj.Timeseries.z = wdenoise(obj.Timeseries.z);
            
        end
        function signal = integrate(obj)
            %
            %
            
            dt = obj.sampleTime.Value/1000;
            int_sig = zeros(3, obj.dataSize);
            int_tim = zeros(1, obj.dataSize);
            int_tim(1,1) = obj.Timeseries.t(1) / 1000;
            for k = 1 : obj.dataSize - 1
                int_tim(1,k+1) = k*dt + int_tim(1,1);
                cur_signal = [obj.Timeseries.x(k); ...
                              obj.Timeseries.y(k); ...
                              obj.Timeseries.z(k)];
                int_sig(:,k+1) = int_sig(:,k) + dt * cur_signal;
            end
            
            signal.x = int_sig(1,:); signal.y = int_sig(2,:);
            signal.z = int_sig(3,:); signal.t = int_tim(1,:);
            
        end
        
        % Show sensor informations
        function summary(obj)
            % This method just prints the summary of the sensor object 
            % information...
            fprintf("<-- Summary of the sensor -> "+ obj.Name +" <-->\n");
            
            fprintf("  (1) TimeSeries collection info:\n");
            fprintf("     |- x size: "+num2str(length(obj.Timeseries.x))+"\n");
            fprintf("     |- y size: "+num2str(length(obj.Timeseries.y))+"\n");
            fprintf("     |- z size: "+num2str(length(obj.Timeseries.z))+"\n");
            fprintf("     |- t size: "+num2str(length(obj.Timeseries.t))+"\n");
            
            fprintf("  (2) Timeseries sampling info:\n");
            fprintf("     |- dT: " + num2str(obj.sampleTime.Value) ...
                + " " + obj.sampleTime.unity + "\n");
            fprintf("     |- dT variation: " + ...
                num2str(obj.sampleTime.Variation) +" "+ ...
                obj.sampleTime.unity+"\n");
            if (obj.digital_data_flag)
                fprintf("     |- This data can be considered digital!\n");
            else
                fprintf("     |- This data cannot be considered digital!\n");
            end
            
            if (obj.resampled_flag)
                fprintf("     |- This data was already resampled!\n");
            end
            fprintf("\n\n");
        end
        function saveState(obj, operation)
           %
           %
           if (isempty(obj.History))
              obj.History = struct(); 
           end
           obj.History.(operation) = obj.Timeseries;
        end
        function plot_serie(obj, varargin)
           % This function just simply plots the sensor time series.
           %
           
           % Default Figure Configs
           cfg.fignum = 1; cfg.lw = 2;
           cfg.xlabel = "Time (s)"; 
           cfg.ylabel = "Time Series Data";
           cfg.legend = ["X","Y","Z"];
           
           if ~isempty(varargin)
               for k = 1 : 2 : length(varargin)
                  cfg.(string(varargin{k})) = varargin{k+1};
               end
           end
            
           data = obj.Timeseries;
           figure(cfg.fignum); hold on;
           plot(data.t/1000, data.x, 'LineWidth',cfg.lw);
           plot(data.t/1000, data.y, 'LineWidth',cfg.lw);
           plot(data.t/1000, data.z, 'LineWidth',cfg.lw);
           legend(cfg.legend); grid on;
           xlabel(cfg.xlabel); ylabel(cfg.ylabel);
           hold off; set(gcf,'color','w');
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

