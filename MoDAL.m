classdef MoDAL
    properties (Constant)
        Version = "1.2.8";
    end

    methods(Static)

        function Help
            fprintf('\nThis class includes a collection of plotting codes often used by MoDAL members.\n\n')
            fprintf('Enter "MoDAL.Install" into the command window to install this class.\n\n')
            fprintf('For help contact Professor Moore at kmoore@unl.edu.\n')

            methods(MoDAL)
        end

        function Install
            p1 = pwd;

            url1 = 'https://drive.google.com/u/0/uc?id=1HOj6OBKCMdXw-BZ5knjRv5dy4e0ED8eq&export=download';
            url2 = 'https://drive.google.com/u/0/uc?id=1aLm6TajZglFrJDN_-XMsgGnAv1UM8TeM&export=download';
            url3 = 'https://drive.google.com/u/0/uc?id=12lB_prP97Tyd1rKCecJcE_7sujuv-axX&export=download';
            url4 = 'https://drive.google.com/u/0/uc?id=1w9ogr_YQAF7pS_suddBKIE_BEesfJGki&export=download';

            filename1 = 'MoDAL.m';
            filename2 = 'emdc_fix.mexw64';
            filename3 = 'emdc_fix.mexmaci64';
            filename4 = 'emdc_fix.m';
            filename5 = 'emdc_fix.mexmaca64';
            source1 = fullfile(p1,filename1);

            destination1 = fullfile([userpath '/'],filename1);
            destination2 = fullfile([userpath '/'],filename2);
            destination3 = fullfile([userpath '/'],filename3);
            destination4 = fullfile([userpath '/'],filename4);
            destination5 = fullfile([userpath '/'],filename5);

            copyfile(source1,destination1,'f')
            if ~isfile(destination2); websave(destination2,url1); end
            if ~isfile(destination3); websave(destination3,url2); end
            if ~isfile(destination4); websave(destination4,url3); end
            if ~isfile(destination5); websave(destination5,url4); end

            destination = fullfile([userpath '/'],'startup.m');

            FID=fopen(destination,'w');
            fprintf(FID,['%%   STARTUPSAV   Startup file\n%%   Change the ' ...
                'name of this file to STARTUP.M. The file\n%%   ' ...
                'is executed when MATLAB starts up, if it exists\n' ...
                '%%   anywhere on the path.  In this example, the\n' ...
                '%%   MAT-file generated during quitting using FINISHSAV\n' ...
                '%%   is loaded into MATLAB during startup.\n\n' ...
                '%%   Copyright 1984-2000 The MathWorks, Inc. \n\n' ...
                'format shortg\nget(0,''Factory'');\n' ...
                'set(0,''defaultfigurecolor'',[1 1 1])\n' ...
                'set(0,''defaultaxesfontsize'',12)\n' ...
                'set(0,''DefaultLineLineWidth'',1);\n' ...
                'MoDAL.Update' ...
                '']);
            fclose(FID);

            fprintf('\nMoDAL version %s successfully installed.\n\n',MoDAL.Version)
        end

        function Update
            URL = 'https://raw.githubusercontent.com/KeeganJMoore/MoDALToolbox/main/MoDAL.m';
            GetRequest = webread(URL);
            R = strfind(GetRequest,'Version');
            P = strfind(GetRequest,';');
            eval(GetRequest(R(1):P(1)));
            if MoDAL.Version ~= Version
                filename1 = 'MoDAL.m';
                destination1 = fullfile([userpath '/'],filename1);
                websave(destination1,URL);
                fprintf('\nMoDAL updated to version %s.\n\n',Version)
            end
        end

        function ExamplePlot
            fprintf(['\nExample for how to use the plotting codes in this object.\nYou can copy-paste the following into a script to play around with the codes.\n\n', ...
                'clc,clear,close all \n', ...
                'time = 0:1e-2:100;\n', ...
                'signal1 = 0.4*sin(2*pi*3.44589*time)+0.138*cos(2*pi*8.48359*time);\n', ...
                'signal2 = 0.3*sin(2*pi*3.44589*time)+0.03578*cos(2*pi*8.48359*time);\n\n', ...
                'T = 1;\n', ...
                'Tb = 3;\n', ...
                'force = sin(pi/T*(time-Tb)).*(heaviside(time-Tb)-heaviside(time-T-Tb));\n', ...
                '%% Plot force\n', ...
                'MoDAL.PlotForce(time,force)\n\n', ...
                'freqMin = 0;\n', ...
                'freqMax = 20;\n\n', ...
                '%% Plot a single sensor\n', ...
                'MoDAL.PlotTSWT(time,signal1,freqMin,freqMax,label=''Disp'')\n', ...
                'MoDAL.PlotTSWTFT(time,signal1,freqMin,freqMax,label=''Disp'')\n\n', ...
                '%% Plot two sensors\n', ...
                'MoDAL.PlotTSWT_Compare(time,signal1,time,signal2,freqMin,freqMax,label=''Disp'')\n', ...
                'MoDAL.PlotTSWTFT_Compare(time,signal1,time,signal2,freqMin,freqMax,label=''Disp'')\n\n'])
            close all
            time = 0:1e-2:100;
            signal1 = 0.4*sin(2*pi*3.44589*time)+0.138*cos(2*pi*8.48359*time);
            signal2 = 0.3*sin(2*pi*3.44589*time)+0.03578*cos(2*pi*8.48359*time);

            T = 1;
            Tb = 3;
            force = sin(pi/T*(time-Tb)).*(heaviside(time-Tb)-heaviside(time-T-Tb));
            force = awgn(force,58);

            % Plot force
            MoDAL.PlotForce(time,force)

            freqMin = 0;
            freqMax = 20;

            % Plot a single sensor
            MoDAL.PlotTSWT(time,signal1,freqMin,freqMax,'label','Disp')
            MoDAL.PlotTSWTFT(time,signal1,freqMin,freqMax,'label','Disp')

            % Plot two sensors
            MoDAL.PlotTSWT_Compare(time,signal1,time,signal2,freqMin,freqMax,'label','Disp')
            MoDAL.PlotTSWTFT_Compare(time,signal1,time,signal2,freqMin,freqMax,'label','Disp')
        end

        function PlotForce(time,force,options)
            %
            % Plots an impact hammer force signal.
            %
            % Required Inputs
            % ---------------------------------------------
            % time - Time vector
            % force - Signal
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % timeStart - lower x-limit for zoomed-in view
            % timeEnd - upper x-limit for zoomed-in view
            % nonDim - 0 sets units to dimensional
            %        - 1 sets units to nondimensional (uses \cdot)
            % forceUnits - character vector that specifies the force units.
            %              nonDim overrides this option.
            % timeUnits - character vector that specifies the time units.
            %             nonDim overrides this option.
            arguments
                time (:,1) double
                force (:,1) double
                options.timeStart double
                options.timeEnd double
                options.nonDim double = 0
                options.forceUnits char = 'N'
                options.timeUnits char = "s"
                options.fontSize = 12;
            end
            tStartSpecified = 1;
            tEndSpecified = 1;
            if options.nonDim == 1
                options.forceUnits = '\cdot';
                options.timeUnits = '\cdot';
            end
            if ~isfield(options,"timeStart")
                idx = findchangepts(force,"MaxNumChanges",2);
                options.timeStart = time(idx(1));
                tStartSpecified = 0;
            end
            if ~isfield(options,"timeEnd")
                idx = findchangepts(force,"MaxNumChanges",2);
                options.timeEnd = time(idx(2));
                tEndSpecified = 0;
            end
            tDiff = (options.timeEnd-options.timeStart)/2;

            if tStartSpecified == 0
                options.timeStart = options.timeStart-tDiff;
            end

            if tEndSpecified == 0
                options.timeEnd = options.timeEnd+tDiff;
            end

            figure
            plot(time,force,'k');
            xlabel(['Time [' options.timeUnits ']'])
            ylabel(['Force [' options.forceUnits ']'])
            xlim([time(1) time(end)])
            ylim(1.1*[min(force) max(force)])
            set(gca,'FontSize',options.fontSize)

            p2 = axes;
            p2.Position = [0.3 0.35 0.55 0.5];
            plot(time,force,'k')
            xlim([options.timeStart options.timeEnd])
            annotation('Arrow','Position',[0.1732,0.2363,0.0786,0.0620]);
            title(sprintf('Max Amplitude = %g %s',round(max(abs(force))),options.forceUnits))
            set(gca,'FontSize',options.fontSize)

        end

        function ProcessRecData(TextFileName,options)
            % This code loads the data stored inside the Run folders. Do not remove the
            % text files from the Run folders.
            %
            % The code assumes that the first channel is the applied force. If this is
            % not true then the code needs to be modified and you should contact Prof.
            % Moore.
            %
            arguments
                TextFileName string
                options.cutOffFreq double = 3;
                options.order double = 3;
                options.endTime double = 1e5;

            end
            A = dir;
            u = 1;
            for ij = 1:length(A)
                if strfind(A(ij).name,'Processed') > 0
                elseif strfind(A(ij).name,'.mat') > 0
                elseif strfind(A(ij).name,'Run') > 0
                    if u == 1
                        FName1 = append(A(ij).name,'/',TextFileName);
                        File1 = fopen(FName1);
                        Qa = textscan(File1,'%s',300);
                        P = char(Qa{1});
                        [~,Wp] = size(P);
                        Zero = '0';
                        for ip = 1:Wp-1
                            Zero = append(Zero,' ');
                        end
                        R = find((sum(P == Zero,2) == Wp) == 1)-1;
                        fclose(File1);
                    end
                    A(ij).name
                    FName1 = append(A(ij).name,'/',TextFileName);
                    File1 = fopen(FName1);
                    Qa = textscan(File1,'%s',R);
                    P = char(Qa{1});
                    [~,Wp] = size(P);
                    Chan = 'Chan';
                    for ip = 1:Wp-4
                        Chan = append(Chan,' ');
                    end
                    N = diff(find((sum(P == Chan,2) == Wp) == 1));
                    C = textscan(File1,'%f',Inf);
                    Aq = C{1,1};
                    fclose(File1);
                    Time = Aq(1:N:end);
                    tB = sum(Time <= options.endTime);
                    Force = detrend(Aq(2:N:end));
                    ForceTemp(:,u) = Force(1:tB);
                    for vp = 3:N
                        AT = detrend(Aq(vp:N:end));
                        AccTemp(:,vp-2,u) = AT(1:tB);
                    end
                    u = u+1;
                end
            end
            Time = Time(1:tB);
            clear Force

            MaxForces = max(abs(ForceTemp));

            Sortu = [(1:u-1)' MaxForces'];
            Sortu = sortrows(Sortu,2);
            for i = 1:u-1
                Force(:,i) = ForceTemp(:,Sortu(i,1));
                Acc(:,:,i) = detrend(AccTemp(:,:,Sortu(i,1)));
            end

            t = Time;
            dt = t(2)-t(1);
            Fs = 1/dt;
            Fnyq = Fs/2;
            Fc = options.cutOffFreq/Fnyq;
            [b,a] = butter(options.order,Fc,'high');

            Data.MaxForce = max(abs(Force))';
            Data.Force = Force;
            Data.Time = Time;
            Data.Acc = Acc;
            Data.Vel = filtfilt(b,a,cumtrapz(t,Data.Acc));
            Data.Disp = filtfilt(b,a,cumtrapz(t,Data.Vel));
            Data.Order = Order;
            Data.CutOffFreq = cutOffFreq;
            Data.TimeEnd = options.endTime;

            save('AllData_Processed.mat','Data')
        end

        % Compute FFT or FRF
        function [f,Fx] = ComputeFFT(time,signal)
            dt = time(2)-time(1);
            L = length(time);
            f = 1/dt*(0:L/2)/L;
            Fx = fft(signal,L)*2/L;
        end

        % Compute WT
        function [freq,modulus] = WaveletTransform(time,signal,minFreq,maxFreq,options)
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Wavelet Transform using Morlet Wavelet
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Inputs
            % ------
            % time - Time vector
            % signal - Signal
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs
            % ---------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            %
            % Outputs
            % -------
            % freq = Frequency vector
            % modulus = Absolute value of wavelet transform of signal
            %
            % Plotting
            % --------
            % This code outputs the modulus of the WT, which represents the energy
            % content at each frequency. Since module is a matrix, there are a couple
            % of ways to plot it. If you are plotting against time, then the best way
            % is to use imagesc as follows:
            %
            % figure
            % imagesc(time,freq,modulus')
            % set(gca,'ydir','nor')
            % colormap(flipud(gray))
            % xlabel('Time [s]')
            % ylabel('Freq. [Hz]')
            %
            % If you are swapping time for amplitude or energy and want to plot on a
            % logarithmic x-scale, then it is best to use pcolor as follows:
            %
            % figure
            % P = pcolor(energy,freq,modulus');
            % P.EdgeColor = 'none';
            % set(gca,'xscale','log')
            % colormap(flipud(gray))
            % xlabel('Energy [J]')
            % ylabel('Freq. [Hz]')
            %
            % Citation
            % --------
            % K.J. Moore, M. Kurt, M. Eriten, D.M. McFarland, L.A. Bergman, A.F. Vakakis,
            % “Wavelet-Bounded Empirical Mode Decomposition for Measured Time Series Analysis,”
            % Mechanical Systems and Signal Processing, 99:14–29, 2018.
            % https://dx.doi.org/10.1016/j.ymssp.2017.06.005

            arguments
                time (:,1) double
                signal (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
            end

            % Transform Parameters
            dt = time(2)-time(1);
            lengthSignal = length(signal);
            df = (maxFreq-minFreq)/options.numFreq; % Frequency Resolution
            freq = minFreq:df:maxFreq; % Prescribe the frequencies of interest
            waveletScale = options.motherWaveletFreq./freq; % Compute the wavelet scales

            % FFT Parameters
            nfourier = 2^nextpow2(lengthSignal);			        % Zero-filling
            npt = nfourier/2;
            Fourierfreq = 1/dt*(0:npt-1)/nfourier; % Frequency vector

            % Compute FFT of xnew
            signalFFT = fft(signal,nfourier);
            signalFFT(npt+1:end) = [];

            % Vectorized Computation of Wavelet Transform of signal
            core2 = bsxfun(@times,conj(bsxfun(@times,(2^0.5)*exp(-0.5* ...
                (2*pi*(bsxfun(@times,Fourierfreq',waveletScale)- ...
                options.motherWaveletFreq)).^2),sqrt(waveletScale))), ...
                signalFFT);

            % Assert Admissibility Condition
            if minFreq == 0
                core2(:,1) = 0;
            end
            result = ifft(core2,nfourier);
            result1 = result(1:lengthSignal,:);
            modulus  = abs(result1);
        end

        % Plot FFT/FRF
        function PlotFFT(time,signal,minFreq,maxFreq,Label,Color,Same)
            % Plots the FFT of the signal.
            % t - Time vector
            % x - Signal
            % minFreq - Minimum frequency shown in the FFT.
            % maxFreq - Maximum frequency shown in the FFT.
            % Label - A character string describing the type of signal (see YLabel).
            % Color - An N x 3 vector for colors.
            % Same - If plotting in a loop, set this to 1 to keep plots on
            %        the same figure.
            arguments
                time double
                signal double
                minFreq double
                maxFreq double
                Label string = '';
                Color double = [0 0 0];
                Same double = 0;
            end

            dt = time(2)-time(1);
            L = length(time);
            f = 1/dt*(0:L/2)/L;
            FFx = 2/L*fft(signal);
            Phase = unwrap(angle(FFx(1:length(f))));
            FFx = abs(FFx(1:length(f)));
            if Same ~= 1
                figure
            end
            subplot(2,1,1)
            semilogy(f,FFx,'Color',Color)
            xlabel('Frequency [Hz]')
            xlim([minFreq maxFreq])
            if strcmp(Label,'Acc')
                ylabel({'Acc.','[m/s^2]'})
            elseif strcmp(Label,'Vel')
                ylabel({'Vel. [m/s]'})
            elseif strcmp(Label,'Disp')
                ylabel({'Disp. [m]'})
            end

            subplot(2,1,2)
            plot(f,Phase,'Color',Color)

            xlabel('Frequency [Hz]')
            ylabel('Angle [rad]')
            xlim([minFreq maxFreq])
        end

        function PlotFRF(time,signal,Force,minFreq,maxFreq,Label,Color,Same)
            % Plots the FRF of the signal.
            % t - Time vector.
            % x - Signal
            % Force -  The applied force vector.
            % minFreq - Minimum frequency shown in the FFT.
            % maxFreq - Maximum frequency shown in the FFT.
            % Label - A character string describing the type of signal (see YLabel).
            % Color - An N x 3 vector for colors.
            % Same - If plotting in a loop, set this to 1 to keep plots on
            %        the same figure.
            arguments
                time (:,1) double
                signal (:,1) double
                Force (:,1) double
                minFreq double
                maxFreq double
                Label string = '';
                Color double = [0 0 0];
                Same double = 0;
            end
            dt = time(2)-time(1);
            L = length(time);
            f = 1/dt*(0:L/2)/L;
            Fx = 2/L*fft(signal,L);
            FF = 2/L*fft(Force,L);
            Phase = Fx(1:length(f))./FF(1:length(f));
            tol = 1e-6;
            Phase(abs(Phase) < tol) = 0;
            Phase = unwrap(angle(Phase));
            FFx = abs(Fx(1:length(f))./FF(1:length(f)));
            if Same ~= 1
                figure
            end
            subplot(2,1,1)
            semilogy(f,FFx,'Color',Color)
            xlabel('Frequency [Hz]')
            xlim([minFreq maxFreq])
            if strcmp(Label,'Acc')
                ylabel({'Accelerance','[(m/s^2)/N]'})
            elseif strcmp(Label,'Vel')
                ylabel({'Mobility','[(m/s)/N]'})
            else
                ylabel({'Receptance','[m/N]'})
            end

            subplot(2,1,2)
            plot(f,Phase,'Color',Color)
            xlabel('Frequency [Hz]')
            ylabel('Angle [rad]')
            xlim([minFreq maxFreq])
        end

        function PlotTSFT(time,signal,minFreq,maxFreq,options)
            % Plots the time series and wavelet transform of a signal.
            %
            % Required Inputs
            % ---------------------------------------------
            % time - Time vector
            % signal - Signal
            % minFreq - Lower limit of the frequency axis for the FFT.
            % maxFreq - Upper limit of the frequency axis for the FFT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % label = label for the plots
            %    'Disp' for displacement data
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            %    'Tension' for tension data
            % timeStart - Sets the minimum xlimit for time plots to this
            %           value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % FontSize - Sets the font size to this value.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first plot.

            arguments
                time (:,1) double
                signal (:,1) double
                minFreq double
                maxFreq double
                options.label string = '';
                options.timeStart double = time(1);
                options.timeEnd double = time(end);
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.force (:,1) double = [];
                options.title string = '';
                options.hideX double = 0;
            end

            figure
            % Plot Time Series
            subplot(2,1,1)
            MoDAL.TSPlot(time,signal,timeStart=options.timeStart,...
                timeEnd=options.timeEnd,label=options.label, ...
                fontSize=options.fontSize,tsLim=options.tsLim)
            title(options.title)

            % Plot FFT
            subplot(2,1,2)
            MoDAL.FTPlot(time,signal,minFreq,maxFreq,force=options.force, ...
                label=options.label,fontSize=options.fontSize)

            if options.hideX; MoDAL.HideX;end
        end

        function PlotTSWT(time,signal,minFreq,maxFreq,options)
            % Plots the time series and wavelet transform of a signal.
            %
            % Required Inputs
            % ---------------------------------------------
            % time - Time vector
            % signal - Signal
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % label = label for the plots
            %    'Disp' for displacement data
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            %    'Tension' for tension data
            % timeStart - Sets the minimum xlimit for time plots to this
            %           value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % FontSize - Sets the font size to this value.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % colorMap - A colormap vector used for coloring the WT
            %            spectrum. Default option is MoDAL.UnbiasedOneMinusPink.
            % wtPower - Raises the WT spectrum to this power. Default
            %           option is 1.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first plot.

            arguments
                time (:,1) double
                signal (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
                options.label string = '';
                options.timeStart double = time(1);
                options.timeEnd double = time(end);
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.colorMap double = MoDAL.UnbiasedOneMinusPink;
                options.wtPower double = 1;
                options.title string = '';
                options.hideX double = 0;
            end

            figure
            % Plot Time Series
            subplot(2,1,1)
            MoDAL.TSPlot(time,signal,timeStart=options.timeStart,...
                timeEnd=options.timeEnd,label=options.label, ...
                fontSize=options.fontSize,tsLim=options.tsLim)
            title(options.title)

            % Compute WT
            [freq,mods] = MoDAL.WaveletSignal(time,signal,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);
            mods = mods/max(mods,[],'All');

            % Plot WT
            subplot(2,1,2)
            MoDAL.WTSpectraPlot(time,freq,mods,options)

            if options.hideX; MoDAL.HideX;end
        end

        function PlotTSWTFT(time,signal,minFreq,maxFreq,options)
            % Plots the time series, wavelet transform, and FFT or FRF of the signal.
            %
            % Required Inputs
            % ---------------------------------------------
            % time - Time vector
            % signal - Signal
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % force – Force signal used to compute FRFs. Defaul is an empty
            %         vector.
            % label = label for the plots
            %    'Disp' for displacement data
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            % timeStart - Sets the minimum xlimit for time plots to this
            %             value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % fontSize - Sets the font size to this value.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % colorMap - A colormap vector used for coloring the WT
            %            spectrum. Default option is MoDAL.UnbiasedOneMinusPink.
            % wtPower - Raises the WT spectrum to this power. Default
            %           option is 1.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first two plots.
            arguments
                time (:,1) double
                signal (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
                options.label string = '';
                options.timeStart double = time(1);
                options.timeEnd double = time(end);
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.force (:,1) double = [];
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.colorMap double = MoDAL.UnbiasedOneMinusPink;
                options.wtPower double = 1;
                options.title string = '';
                options.hideX double = 0;
            end

            figure
            % Plot Time Series
            subplot(3,1,1)
            MoDAL.TSPlot(time,signal,timeStart=options.timeStart,...
                timeEnd=options.timeEnd,label=options.label, ...
                fontSize=options.fontSize,tsLim=options.tsLim)
            title(options.title)

            % Compute WT
            [freq,mods] = MoDAL.WaveletSignal(time,signal,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);
            mods = mods/max(mods,[],'All');

            % Plot WT
            subplot(3,1,2)
            MoDAL.WTSpectraPlot(time,freq,mods,options)

            % Plot FFT/FRF
            subplot(3,1,3)
            MoDAL.FTPlot(time,signal,minFreq,maxFreq,force=options.force, ...
                label=options.label,fontSize=options.fontSize)

            if options.hideX; MoDAL.HideX;end
        end

        function PlotTSMWT(time,signal,minFreq,maxFreq,options)
            % Plots the time series and the MWT of the signal.
            %
            % Required Inputs
            % ---------------------------------------------
            % time - Time vector
            % signal - Signal
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % label = label for the plots
            %    'Disp' for displacement data
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            % timeStart - Sets the minimum xlimit for time plots to this
            %             value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % fontSize - Sets the font size to this value.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first plot.
            arguments
                time (:,1) double
                signal (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
                options.label string = '';
                options.timeStart double = time(1);
                options.timeEnd double = time(end);
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.title string = '';
                options.hideX double = 0;
            end

            figure
            % Plot Time Series
            subplot(2,1,1)
            MoDAL.TSPlot(time,signal,timeStart=options.timeStart,...
                timeEnd=options.timeEnd,label=options.label)
            title(options.title)

            % Compute WT
            [freq,mods] = MoDAL.WaveletSignal(time,signal,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);
            mods = mods/max(mods,[],'All');

            % Plot MWT
            subplot(2,1,2)
            MoDAL.MWTPlot(freq,max(mods),minFreq,maxFreq, ...
                fontSize=options.fontSize,label=options.label)

            if options.hideX; MoDAL.HideX;end
        end

        function PlotTSWTMWT(time,signal,minFreq,maxFreq,options)
            % Plots the time series, wavelet transform, and the MWT of the signal.
            %
            % Required Inputs
            % ---------------------------------------------
            % time - Time vector
            % signal - Signal
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % label = label for the plots
            %    'Disp' for displacement data
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            % timeStart - Sets the minimum xlimit for time plots to this
            %             value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % fontSize - Sets the font size to this value.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % colorMap - A colormap vector used for coloring the WT
            %            spectrum. Default option is MoDAL.UnbiasedOneMinusPink.
            % wtPower - Raises the WT spectrum to this power. Default
            %           option is 1.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first two plots.
            arguments
                time (:,1) double
                signal (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
                options.label string = '';
                options.timeStart double = time(1);
                options.timeEnd double = time(end);
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.colorMap double = MoDAL.UnbiasedOneMinusPink;
                options.wtPower double = 1;
                options.title string = '';
                options.hideX double = 0;
            end

            figure
            % Plot Time Series
            subplot(3,1,1)
            MoDAL.TSPlot(time,signal,timeStart=options.timeStart,...
                timeEnd=options.timeEnd,label=options.label)
            title(options.title)

            % Compute the WT
            [freq,mods] = MoDAL.WaveletSignal(time,signal,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);
            mods = mods/max(mods,[],'All');

            % Plot the WT
            subplot(3,1,2)
            MoDAL.WTSpectraPlot(time,freq,mods,options)

            % Plot the MWT
            subplot(3,1,3)
            MoDAL.MWTPlot(freq,max(mods),minFreq,maxFreq, ...
                fontSize=options.fontSize,label=options.label)

            if options.hideX; MoDAL.HideX;end
        end

        function PlotTSWT_Compare(time1,signal1,time2,signal2,minFreq,maxFreq,options)
            % Plots the time series and wavelet transform of two signals.
            %
            % Required Inputs
            % ---------------------------------------------
            % time1 = time vector for first time series
            % signal1 = first time series
            % time2 = time vector for second time series
            % signal2 = second time series
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % label = label for the plots
            %    'Disp' for displacement data
            %    'DispND' for non-dimensional displacement
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            %    'Tension' for tension data
            % timeStart - Sets the minimum xlimit for time plots to this
            %             value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single 
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % fontSize - Sets the font size to this value.
            % firstColor - Color string for the first signal.
            % secondColor - Color string for the second signal.
            % firstLineStyle – Line style for the first signal.
            % secondLineStyle – Line style for the first signal.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % colorMap - A colormap vector used for coloring the WT
            %            spectrum. Default option is MoDAL.UnbiasedOneMinusPink.
            % legends - 2x1 cell containing two strings that describe
            %           signal 1 and signal 2, respectively.
            % wtPower - Raises the WT spectrum to this power. Default
            %           option is 1.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first two plots.
            arguments
                time1 (:,1) double
                signal1 (:,1) double
                time2 (:,1) double
                signal2 (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
                options.label string = '';
                options.timeStart double = min(time1(1),time2(1));
                options.timeEnd double = max(time1(end),time2(end));
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.firstColor string = 'k';
                options.secondColor string = 'r';
                options.firstLineStyle string = '-';
                options.secondLineStyle string = '-';
                options.force double = [];
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.colorMap double = MoDAL.UnbiasedOneMinusPink;
                options.legends cell = {'Signal 1','Signal 2'};
                options.legendsLoc string = 'northeast';
                options.wtPower double = 1;
                options.title string = '';
                options.hideX double = 0;
            end

            figure
            % Plot Time Series
            subplot(3,1,1);
            MoDAL.TSPlot(time1,signal1,color=options.firstColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.firstLineStyle)
            hold on
            MoDAL.TSPlot(time2,signal2,color=options.secondColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.secondLineStyle)
            legend(options.legends,'Location',options.legendsLoc,'orientation','horizontal')
            title(options.title)

            % Compute Wavelet Transforms
            [freq1,mods1] = MoDAL.WaveletSignal(time1,signal1,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);
            [freq2,mods2] = MoDAL.WaveletSignal(time2,signal2,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);

            % Normalize Wavelets
            maxnorm = max([max(max(mods1)) max(max(mods2))]);
            mods1 = mods1/maxnorm;
            mods2 = mods2/maxnorm;

            % Plot Wavelets
            subplot(3,1,2);
            MoDAL.WTSpectraPlot(time1,freq1,mods1,options)
            clim([0 1])
            title(options.legends{1})

            subplot(3,1,3);
            MoDAL.WTSpectraPlot(time2,freq2,mods2,options)
            clim([0 1])
            title(options.legends{2})
            if options.hideX; MoDAL.HideX;end
        end


        function PlotTSWTFT_Compare(time1,signal1,time2,signal2,minFreq,maxFreq,options)
            % Plots the time series, wavelet transform, and FFT/FRF of two signals.
            %
            % Required Inputs
            % ---------------------------------------------
            % time1 = time vector for first time series
            % signal1 = first time series
            % time2 = time vector for second time series
            % signal2 = second time series
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % force – Force signal used to compute FRFs. Defaul is an empty
            %         vector.
            % label = label for the plots
            %    'Disp' for displacement data
            %    'DispND' for non-dimensional displacement
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            %    'Tension' for tension data
            % timeStart - Sets the minimum xlimit for time plots to this
            %             value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % fontSize - Sets the font size to this value.
            % firstColor - Color string for the first signal.
            % secondColor - Color string for the second signal.
            % firstLineStyle – Line style for the first signal.
            % secondLineStyle – Line style for the first signal.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % colorMap - A colormap vector used for coloring the WT
            %            spectrum. Default option is MoDAL.UnbiasedOneMinusPink.
            % legends - 2x1 cell containing two strings that describe
            %           signal 1 and signal 2, respectively.
            % wtPower - Raises the WT spectrum to this power. Default
            %           option is 1.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first two plots.
            arguments
                time1 (:,1) double
                signal1 (:,1) double
                time2 (:,1) double
                signal2 (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
                options.label string = '';
                options.timeStart double = min(time1(1),time2(1));
                options.timeEnd double = max(time1(end),time2(end));
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.firstColor string = 'k';
                options.secondColor string = 'r';
                options.firstLineStyle string = '-';
                options.secondLineStyle string = '-';
                options.force double = [];
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.colorMap double = MoDAL.UnbiasedOneMinusPink;
                options.legends cell = {'Signal 1','Signal 2'};
                options.legendsLoc string = 'northeast';
                options.wtPower double = 1;
                options.title string = '';
                options.hideX double = 0;
            end
            if size(options.force,1) < size(options.force,2)
                options.force = options.force';
            end
            if size(options.force,2) > 1
                F1 = options.force(:,1);
                F2 = options.force(:,2);
            else
                F1 = options.force;
                F2 = options.force;
            end

            figure
            subplot(3,1,1)
            MoDAL.TSPlot(time1,signal1,color=options.firstColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.firstLineStyle)
            hold on
            MoDAL.TSPlot(time2,signal2,color=options.secondColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.secondLineStyle)
            legend(options.legends,'Location',options.legendsLoc,'orientation','horizontal')
            title(options.title)

            % Compute Wavelet Transforms
            [freq1,mods1] = MoDAL.WaveletSignal(time1,signal1,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);
            [freq2,mods2] = MoDAL.WaveletSignal(time2,signal2,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);

            % Normalize Wavelets
            maxnorm = max([max(max(mods1)) max(max(mods2))]);
            mods1 = mods1/maxnorm;
            mods2 = mods2/maxnorm;

            % Plot Wavelets
            subplot(3,2,3)
            MoDAL.WTSpectraPlot(time1,freq1,mods1,options)
            title(options.legends{1})
            clim([0 1])

            subplot(3,2,4)
            MoDAL.WTSpectraPlot(time2,freq2,mods2,options)
            title(options.legends{2})
            clim([0 1])

            S4 = subplot(3,1,3);
            MoDAL.FTPlot(time1,signal1,minFreq,maxFreq,force=F1,color=options.firstColor, ...
                label=options.label,fontSize=options.fontSize, ...
                linestyle=options.firstLineStyle)
            hold on
            MoDAL.FTPlot(time2,signal2,minFreq,maxFreq,force=F2,color=options.secondColor, ...
                label=options.label,fontSize=options.fontSize,linestyle=options.secondLineStyle)
            legend(options.legends,'Location',options.legendsLoc,'orientation','horizontal')
            drawnow;
            S4.Position(2) = 0.12464;
            if options.hideX; MoDAL.HideX;end
        end

        function PlotTSMWT_Compare(time1,signal1,time2,signal2,minFreq,maxFreq,options)
            % Plots the time series and MWT of two signals.
            %
            % Required Inputs
            % ---------------------------------------------
            % time1 = time vector for first time series
            % signal1 = first time series
            % time2 = time vector for second time series
            % signal2 = second time series
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % label = label for the plots
            %    'Disp' for displacement data
            %    'DispND' for non-dimensional displacement
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            %    'Tension' for tension data
            % timeStart - Sets the minimum xlimit for time plots to this
            %             value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % fontSize - Sets the font size to this value.
            % firstColor - Color string for the first signal.
            % secondColor - Color string for the second signal.
            % firstLineStyle – Line style for the first signal.
            % secondLineStyle – Line style for the first signal.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % legends - 2x1 cell containing two strings that describe
            %           signal 1 and signal 2, respectively.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first two plots.
            arguments
                time1 (:,1) double
                signal1 (:,1) double
                time2 (:,1) double
                signal2 (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
                options.label string = '';
                options.timeStart double = min(time1(1),time2(1));
                options.timeEnd double = max(time1(end),time2(end));
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.firstColor string = 'k';
                options.secondColor string = 'r';
                options.firstLineStyle string = '-';
                options.secondLineStyle string = '-';
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.legends cell = {'Signal 1','Signal 2'};
                options.legendsLoc string = 'northeast';
                options.title string = '';
                options.hideX double = 0;
            end

            figure

            % Plot the Time Series
            subplot(2,1,1)
            MoDAL.TSPlot(time1,signal1,color=options.firstColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.firstLineStyle)
            hold on
            MoDAL.TSPlot(time2,signal2,color=options.secondColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.secondLineStyle)
            legend(options.legends,'Location',options.legendsLoc,'orientation','horizontal')
            title(options.title)

            % Compute Wavelet Transforms
            [freq1,mods1] = MoDAL.WaveletSignal(time1,signal1,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);
            [freq2,mods2] = MoDAL.WaveletSignal(time2,signal2,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);

            % Normalize Wavelets
            maxnorm = max([max(max(mods1)) max(max(mods2))]);
            mods1 = mods1/maxnorm;
            mods2 = mods2/maxnorm;

            % Plot MWT
            subplot(2,1,2);
            MoDAL.MWTPlot(freq1,max(mods1),minFreq,maxFreq, ...
                fontSize=options.fontSize,label=options.label)
            hold on
            MoDAL.MWTPlot(freq2,max(mods2),minFreq,maxFreq,color='c', ...
                fontSize=options.fontSize,label=options.label)
            legend(options.legends,'Location',options.legendsLoc, ...
                'orientation','horizontal')

            if options.hideX; MoDAL.HideX;end
            drawnow;
        end

        function PlotTSWTMWT_Compare(time1,signal1,time2,signal2,minFreq,maxFreq,options)
            % Plots the time series, wavelet transform, and MWT of two signals.
            %
            % Required Inputs
            % ---------------------------------------------
            % time1 = time vector for first time series
            % signal1 = first time series
            % time2 = time vector for second time series
            % signal2 = second time series
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % label = label for the plots
            %    'Disp' for displacement data
            %    'DispND' for non-dimensional displacement
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            %    'Tension' for tension data
            % timeStart - Sets the minimum xlimit for time plots to this
            %             value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % fontSize - Sets the font size to this value.
            % firstColor - Color string for the first signal.
            % secondColor - Color string for the second signal.
            % firstLineStyle – Line style for the first signal.
            % secondLineStyle – Line style for the first signal.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % colorMap - A colormap vector used for coloring the WT
            %            spectrum. Default option is MoDAL.UnbiasedOneMinusPink.
            % legends - 2x1 cell containing two strings that describe
            %           signal 1 and signal 2, respectively.
            % wtPower - Raises the WT spectrum to this power. Default
            %           option is 1.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first two plots.
            arguments
                time1 (:,1) double
                signal1 (:,1) double
                time2 (:,1) double
                signal2 (:,1) double
                minFreq double
                maxFreq double
                options.numFreq double = 100;
                options.motherWaveletFreq double = 2;
                options.label string = '';
                options.timeStart double = min(time1(1),time2(1));
                options.timeEnd double = max(time1(end),time2(end));
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.firstColor string = 'k';
                options.secondColor string = 'r';
                options.firstLineStyle string = '-';
                options.secondLineStyle string = '-';
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.colorMap double = MoDAL.UnbiasedOneMinusPink;
                options.legends cell = {'Signal 1','Signal 2'};
                options.legendsLoc string = 'northeast';
                options.wtPower double = 1;
                options.title string = '';
                options.hideX double = 0;
            end

            figure
            subplot(3,1,1)
            MoDAL.TSPlot(time1,signal1,color=options.firstColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.firstLineStyle)
            hold on
            MoDAL.TSPlot(time2,signal2,color=options.secondColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.secondLineStyle)
            legend(options.legends,'Location',options.legendsLoc,'orientation','horizontal')
            title(options.title)

            % Compute Wavelet Transforms
            [freq1,mods1] = MoDAL.WaveletSignal(time1,signal1,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);
            [freq2,mods2] = MoDAL.WaveletSignal(time2,signal2,minFreq,maxFreq, ...
                options.numFreq,options.motherWaveletFreq,options.mirrori,options.mirrorf);

            % Normalize Wavelets
            maxnorm = max([max(max(mods1)) max(max(mods2))]);
            mods1 = mods1/maxnorm;
            mods2 = mods2/maxnorm;

            % Plot Wavelets
            subplot(3,2,3)
            MoDAL.WTSpectraPlot(time1,freq1,mods1,options)
            title(options.legends{1})

            subplot(3,2,4)
            MoDAL.WTSpectraPlot(time2,freq2,mods2,options)
            title(options.legends{2})

            % Plot MWT
            S4 = subplot(3,1,3);
            MoDAL.MWTPlot(freq1,max(mods1),minFreq,maxFreq, ...
                fontSize=options.fontSize,label=options.label,color=options.firstColor)
            hold on
            MoDAL.MWTPlot(freq2,max(mods2),minFreq,maxFreq,color=options.secondColor, ...
                fontSize=options.fontSize,label=options.label)
            legend(options.legends,'Location',options.legendsLoc, ...
                'orientation','horizontal')
            drawnow;
            S4.Position(2) = 0.12464;
            if options.hideX; MoDAL.HideX;end
        end


        function PlotTSFT_Compare(time1,signal1,time2,signal2,minFreq,maxFreq,options)
            % Plots the time series, wavelet transform, and FFT/FRF of two signals.
            %
            % Inputs
            % ---------------------------------------------
            % time1 = time vector for first time series
            % signal1 = first time series
            % time2 = time vector for second time series
            % signal2 = second time series
            % minFreq - Minimum frequency computed in the WT.
            % maxFreq - Maximum frequency computed in the WT.
            %
            % Options
            % -------
            % numFreq - The number of frequency points used in the WT
            %           (default value = 100).
            % motherWaveletFreq - Frequency of the mother wavelet used in
            %           the WT (default value = 2).
            % label = label for the plots
            %    'Disp' for displacement data
            %    'DispND' for non-dimensional displacement
            %    'Vel' for velocity data
            %    'Acc' for acceleration data
            %    'Force' for force data
            %    'Modal' for modal displacement data
            %    'Tension' for tension data
            % timeStart - Sets the minimum xlimit for time plots to this
            %             value.
            % timeEnd - Sets the maximum xlimit for time plots to this
            %           value.
            % tsLim - Sets the ylimit to [-tsLim, tsLim] if only a single
            %         value is provided or to [tsLim(1) tsLim(2)] if two are
            %         provided.
            % fontSize - Sets the font size to this value.
            % firstColor - Color string for the first signal.
            % secondColor - Color string for the second signal.
            % firstLineStyle – Line style for the first signal.
            % secondLineStyle – Line style for the first signal.
            % mirrori - Mirroring applied to the beginning of the signal.
            %           Even mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % mirrorf - Mirroring applied to the end of the signal. Even
            %           mirroring is the default option. Input should
            %           be a string with either 'e', 'o', or 'none'.
            % colorMap - A colormap vector used for coloring the WT
            %            spectrum. Default option is MoDAL.UnbiasedOneMinusPink.
            % legends - 2x1 cell containing two strings that describe
            %           signal 1 and signal 2, respectively.
            % wtPower - Raises the WT spectrum to this power. Default
            %           option is 1.
            % title - Adds title to the time series plot. Default is no
            %         title.
            % hideX - Hides the xticklabels and xlabel for the first two plots.
            arguments
                time1 (:,1) double
                signal1 (:,1) double
                time2 (:,1) double
                signal2 (:,1) double
                minFreq double
                maxFreq double
                options.label string = '';
                options.timeStart double = min(time1(1),time2(1));
                options.timeEnd double = max(time1(end),time2(end));
                options.tsLim double = nan;
                options.fontSize double = 12;
                options.firstColor string = 'k';
                options.secondColor string = 'r';
                options.firstLineStyle string = '-';
                options.secondLineStyle string = '-';
                options.force double = [];
                options.legends cell = {'Signal 1','Signal 2'};
                options.legendsLoc string = 'northeast';
                options.title string = '';
                options.hideX double = 0;
            end
            if size(options.force,1) < size(options.force,2)
                options.force = options.force';
            end
            if size(options.force,2) > 1
                F1 = options.force(:,1);
                F2 = options.force(:,2);
            else
                F1 = options.force;
                F2 = options.force;
            end

            figure
            subplot(2,1,1)
            MoDAL.TSPlot(time1,signal1,color=options.firstColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.firstLineStyle)
            hold on
            MoDAL.TSPlot(time2,signal2,color=options.secondColor,fontSize=options.fontSize, ...
                label=options.label,timeStart=options.timeStart,timeEnd=options.timeEnd, ...
                tsLim=options.tsLim,linestyle=options.secondLineStyle)
            legend(options.legends,'Location',options.legendsLoc,'orientation','horizontal')
            title(options.title)

            subplot(2,1,2);
            MoDAL.FTPlot(time1,signal1,minFreq,maxFreq,force=F1,color=options.firstColor, ...
                label=options.label,fontSize=options.fontSize,linestyle=options.firstLineStyle)
            hold on
            MoDAL.FTPlot(time2,signal2,minFreq,maxFreq,force=F2,color=options.secondColor, ...
                label=options.label,fontSize=options.fontSize,linestyle=options.secondLineStyle)
            legend(options.legends,'Location',options.legendsLoc,'orientation','horizontal')
            drawnow;
            if options.hideX; MoDAL.HideX;end
        end

        function out = UnbiasedOneMinusPink
            out = MoDAL.UnbiasedOneMinusPinkMat;
        end

        % Plot FE Models
        function PlotMesh(coordinates,nodes,options)
            %--------------------------------------------------------------------------
            % Plots a finite element mesh
            %
            % Required Inputs
            % ---------------------------------------------
            % coordinates - The nodal coordinates of the mesh
            %               -----> coordinates = [X Y Z]
            % nodes - The nodal connectivity of the elements
            %         -----> nodes = [node1 node2......]
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % showNodeNumber - display node numbers:
            %                  0 (default) - do not display
            %                  1           - display
            % showFill - include fill in the elements
            %            0 (default) - do not include
            %            1 - include
            % coordzero -
            %
            %
            % Original version coded by:
            %    Siva Srinivas Kolukula, PhD
            %    Indian Tsunami Early Warning Centre (ITEWC)
            %    Advisory Services and Satellite Oceanography Group (ASG)
            %    Indian National Centre for Ocean Information Services (INCOIS)
            %    Hyderabad, INDIA
            %    E-mail: allwayzitzme@gmail.com
            %    web-link: https://sites.google.com/site/kolukulasivasrinivas/
            %
            %    Version 1: 28 August 2011
            %    Version 2: 16 September 2016
            %
            % Modified by Keegan Moore on 25 October 2023

            arguments
                coordinates double
                nodes double
                options.coordzero double = [];
                options.showNodeNumber double = 0;
                options.showFill double = 0;
                options.showElementNumber double = 0;
                options.edgeColor = 'k';
            end

            dimension = size(coordinates,2);  % Dimension of the mesh
            nel = length(nodes) ;                  % number of elements
            nnode = length(coordinates) ;          % total number of nodes in system
            nnel = size(nodes,2);                % number of nodes per element

            % Initialization of the required matrices
            X = zeros(nnel,nel);
            Y = X;
            Z = X;
            Z2 = X;

            if dimension == 3   % For 3D plots
                if nnel==4 % surface in 3D
                    for iel=1:nel
                        nd = nodes(iel,:) ;
                        X(:,iel)=coordinates(nd,1);    % extract x value of the node
                        Y(:,iel)=coordinates(nd,2);    % extract y value of the node
                        Z(:,iel)=coordinates(nd,3) ;   % extract z value of the node
                        if ~isempty(options.coordzero)
                            Z2(:,iel)=options.coordzero(nd);
                        end
                    end
                    if isempty(options.coordzero)
                        Z2 = Z;
                    end
                    if options.showFill == 0
                        patch(X,Y,Z,Z,'facecolor','none','edgecolor',options.edgeColor)
                    else
                        patch(X,Y,Z,Z2,'edgecolor',options.edgeColor)
                    end

                    % display Node numbers and Element numbers
                    if options.showNodeNumber ~= 0
                        for i = 1:length(coordinates(:,1))
                            text(coordinates(i,1),coordinates(i,2),coordinates(i,3),int2str(i),....
                                'fontsize',10,'color','r');
                        end
                    end
                    if options.showElementNumber ~= 0
                        for i = 1:nel
                            text(sum(X(:,i))/4,sum(Y(:,i))/4,sum(Z(:,i))/4,int2str(i),.....
                                'fontsize',10,'color','g') ;
                        end
                    end

                elseif nnel==8  % solid in 3D
                    fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
                    XYZ = cell(1,nel) ;
                    for e=1:nel
                        nd=nodes(e,:);
                        X(:,e) = coordinates(nd,1) ;
                        Y(:,e) = coordinates(nd,2) ;
                        Z(:,e) = coordinates(nd,3) ;
                        XYZ{e} = [X(:,e)  Y(:,e) Z(:,e)] ;
                    end
                    % Plot FEM mesh

                    cellfun(@patch,repmat({'Vertices'},1,nel),XYZ,.......
                        repmat({'Faces'},1,nel),repmat({fm},1,nel),......
                        repmat({'FaceColor'},1,nel),repmat({'r'},1,nel));

                    set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;

                    % display Node numbers and Element numbers
                    if options.showElementNumber ~= 0
                        for i = 1:nel
                            text(X(:,i),Y(:,i),Z(:,i),int2str(nd(i)),....
                                'fontsize',8,'color','k');
                            text(sum(X(:,i))/8,sum(Y(:,i))/8,sum(Z(:,i))/8,int2str(i),.....
                                'fontsize',10,'color','r') ;
                        end
                    end
                end

            elseif dimension == 2           % For 2D plots
                for iel=1:nel
                    nd = nodes(iel,:) ;
                    X(:,iel)=coordinates(nd,1);    % extract x value of the node
                    Y(:,iel)=coordinates(nd,2);    % extract y value of the node
                end

                % Plotting the FEM mesh, diaplay Node numbers and Element numbers
                fill(X,Y,'w')

                if options.showElementNumber ~= 0
                    k = 1:nnode ;
                    nd = k' ;
                    for i = 1:nel
                        text(X(:,i),Y(:,i),int2str(nd(i)),'fontsize',8,'color','k');
                        text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
                    end
                end
            end
            view(3)

        end

        %%%%%%% Decomposition Routines %%%%%%%%
        function Mode = IFD(time,signal,lowerFreqs,upperFreqs)
            dt = time(2)-time(1);

            %FFT Parameters
            ntemps = length(time);
            power2 = nextpow2(ntemps);
            L = 2^power2;
            f = 1/dt*(0:L/2)/L;

            % Compute FFT of x
            xFFT = fft(signal,L);

            Mode = zeros(length(time),length(lowerFreqs));
            for i = 1:length(lowerFreqs)
                FA = sum(f <= lowerFreqs(i));
                FB = sum(f <= upperFreqs(i));

                F = zeros(L/2,1);
                F(FA:FB) = 1;
                F = [F;flipud(F)];
                F = F(1:L);
                fft_mode = F.*xFFT;
                ModeS = ifft(fft_mode,L,'symmetric');
                Mode(:,i) = ModeS(1:length(time));
            end
        end

        function Mode = IWD(time,signal,lowerFreqs,upperFreqs,options)
            arguments
                time double
                signal double
                lowerFreqs double
                upperFreqs double
                options.motherWaveletFreq = 16;
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.chp1 = 0.2;
                options.chp2 = 0.2;
            end

            [x_mirror,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd] = ...
                MoDAL.MirrorSignal(time,signal,options.mirrori, ...
                options.mirrorf,options.chp1,options.chp2);
            
            signal = x_mirror;
        
            Fo = options.motherWaveletFreq;
            dt = time(2)-time(1);
            Fs = 1/dt;

            l_new = 2^nextpow2(size(signal,1));
            nfourier = l_new;
            npt = nfourier/2;
            freq = Fs*(0:nfourier-1)/nfourier;
            FREQ = Fs*([0:npt-1 npt:-1:1])/nfourier;
            interval_freq = FREQ;
            a = Fo./interval_freq;
            fft_MW = conj(bsxfun(@times,pi^(1/4)*(2^0.5)*(exp(-0.5*(2*pi*...
                (bsxfun(@times,FREQ',a)-Fo)).^2)-1*exp(-0.5*(2*pi^2*...
                (bsxfun(@times,FREQ',a).^2+Fo.^2)))),sqrt(a)));

            FA = sum(freq <= lowerFreqs);
            FB = sum(freq <= upperFreqs);

            tff = fft(signal,nfourier);
            noyau2 = bsxfun(@times,fft_MW,tff);
            CMod = ifft(noyau2,nfourier);

            ModeWT = 0*CMod;
            ModeWT(:,FA:FB) = CMod(:,FA:FB);

            fft_Mode1 = fft(ModeWT,nfourier);
            fft_mode_2D = fft_Mode1./fft_MW;
            fft_mode = diag(fft_mode_2D);
            fft_mode(isnan(fft_mode)) = 0;
            Mode = ifft(fft_mode,nfourier,'symmetric');
            if NoMirrorIni == 0 && NoMirrorEnd == 0
                Mode = Mode(L_chp1:end-L_chp2+1);
            elseif NoMirrorIni == 0 && NoMirrorEnd == 1
                Mode = Mode(L_chp1:end);
            elseif NoMirrorIni == 1 && NoMirrorEnd == 0
                Mode = Mode(1:end-L_chp2+1);
            end
            Mode = Mode(1:length(time));
        end

        function Mode = IWD_Mult(time,signal,lowerFreqs,upperFreqs)
            Fo = 2;
            dt = time(2)-time(1);
            Fs = 1/dt;

            l_new = 2^nextpow2(size(signal,1));
            nfourier = l_new;
            npt = nfourier/2;
            freq = Fs*(0:(nfourier-1))/nfourier;
            FREQ = Fs*([0:npt-1 npt:-1:1])/nfourier;
            interval_freq = FREQ;
            a = Fo./interval_freq;
            fft_MW = conj(bsxfun(@times,pi^(1/4)*(2^0.5)*(exp(-0.5* ...
                (2*pi*(bsxfun(@times,FREQ',a)-Fo)).^2)-1*exp(-0.5* ...
                (2*pi^2*(bsxfun(@times,FREQ',a).^2+Fo.^2)))),sqrt(a)));

            tff = fft(signal,nfourier);
            noyau2 = bsxfun(@times,fft_MW,tff);
            CMod = ifft(noyau2,nfourier);
            Mode = zeros(length(time),length(lowerFreqs));
            for i = 1:length(lowerFreqs)
                FA = sum(freq <= lowerFreqs(i));
                FB = sum(freq <= upperFreqs(i));
                ModeWT = 0*CMod;
                ModeWT(:,FA:FB) = CMod(:,FA:FB);

                fft_Mode1 = fft(ModeWT,nfourier);
                fft_mode_2D = fft_Mode1./fft_MW;
                fft_mode = diag(fft_mode_2D);
                fft_mode(isnan(fft_mode)) = 0;
                ModeS = ifft(fft_mode,nfourier,'symmetric');
                Mode(:,i) = ModeS(1:length(time));
            end
        end

        function YLabel(Label,NF)
            % Adds a label to the yaxis of the current axis. Accepts custom labels or a predefined,
            % case insensitive label. The predefined labels are:
            %
            % disp          =>   Disp. [m]
            % dispmm        =>   Disp. [mm]
            % dispin        =>   Disp. [in]
            % dispft        =>   Disp. [ft]
            % dispnd        =>   Disp. [\cdot] (non-dimensional)
            % vel           =>   Vel. [m/s]
            % velmm         =>   Vel. [mm/s]
            % velin         =>   Vel. [in/s]
            % velft         =>   Vel. [ft/s]
            % velnd         =>   Vel. [\cdot] (non-dimensional)
            % accmm         =>   Acc. [mm/s^2]
            % acc           =>   Acc. [m/s^2]
            % accin         =>   Acc. [in/s^2]
            % accft         =>   Acc. [ft/s^2]
            % accnd         =>   Acc. [\cdot] (non-dimensional)
            % force         =>   Force [N]
            % forcekn       =>   Force [kN]
            % forcelb       =>   Force [lb]
            % forcekip      =>   Force [kip]
            % forcend       =>   Force [\cdot] (non-dimensional)
            % tension       =>   Tension [N]
            % tensionkn     =>   Tension [kN]
            % tensionlb     =>   Tension [lb]
            % tensionkip    =>   Tension [kip]
            % tensionnd     =>   Tension [\cdot] (non-dimensional)
            % modal         =>   Amplitude [\cdot]
            % strain        =>   Strain [m/m]
            % microstrain   =>   Strain [\mum/m]
            % strainin      =>   Strain [in/in]
            % strainft      =>   Strain [ft/ft]
            % imf           =>   Amplitude [\cdot]
            % imfdm         =>   Amplitude [m]
            % imfdmm        =>   Amplitude [mm]
            % imfdin        =>   Amplitude [in]
            % imfdft        =>   Amplitude [ft]
            % imfvm         =>   Amplitude [m/s]
            % imfvmm        =>   Amplitude [mm/s]
            % imfvin        =>   Amplitude [in/s]
            % imfvft        =>   Amplitude [ft/s]
            % imfam         =>   Amplitude [m/s^2]
            % imfamm        =>   Amplitude [mm/s^2]
            % imfain        =>   Amplitude [in/s^2]
            % imfaft        =>   Amplitude [ft/s^2]
            % imff          =>   Amplitude [N]
            % imffn         =>   Amplitude [kN]
            % imfflb        =>   Amplitude [lb]
            % imffkip       =>   Amplitude [kip]
            % phase         =>   Phase [rad]
            % phasevel      =>   Phase Vel. [rad/s]
            
            if nargin == 1
                switch lower(Label)
                    case 'disp'
                        ylabel('Disp. [m]')
                    case 'dispin'
                        ylabel('Disp. [in]')
                    case 'dispft'
                        ylabel('Disp. [in]')
                    case 'dispmm'
                        ylabel('Disp. [mm]')
                    case 'dispnd'
                        ylabel('Disp. [\cdot]')
                    case 'vel'
                        ylabel('Vel. [m/s]')
                    case 'velmm'
                        ylabel('Vel. [mm/s]')
                    case 'velin'
                        ylabel('Vel. [in/s]')
                    case 'velft'
                        ylabel('Vel. [ft/s]')
                    case 'velnd'
                        ylabel('Vel. [\cdot]')
                    case 'acc'
                        ylabel('Acc. [m/s^2]')
                    case 'accmm'
                        ylabel('Acc. [mm/s^2]')
                    case 'accin'
                        ylabel('Acc. [in/s^2]')
                    case 'accft'
                        ylabel('Acc. [ft/s^2]')
                    case 'accnd'
                        ylabel('Acc. [\cdot]')
                    case 'force'
                        ylabel('Force [N]')
                    case 'forcekn'
                        ylabel('Force [kN]')
                    case 'forcelb'
                        ylabel('Force [lb]')
                    case 'forcekip'
                        ylabel('Force [kip]')
                    case 'forcend'
                        ylabel('Force [\cdot]')
                    case 'tension'
                        ylabel('Tension [N]')
                    case 'tensionkn'
                        ylabel('Tension [kN]')
                    case 'tensionlb'
                        ylabel('Tension [lb]')
                    case 'tensionkip'
                        ylabel('Tension [kip]')
                    case 'tensionnd'
                        ylabel('Tension [\cdot]')
                    case 'modal'
                        ylabel('Amplitude [\cdot]')
                    case 'strain'
                        ylabel('Strain [m/m]')
                    case 'microstrain'
                        ylabel('Strain [\mum/m]')
                    case 'strainin'
                        ylabel('Strain [in/in]')
                    case 'strainft'
                        ylabel('Strain [ft/ft]')
                    case 'imf'
                        ylabel('Amplitude [\cdot]')
                    case 'imfdm'
                        ylabel('Amplitude [m]')
                    case 'imfdmm'
                        ylabel('Amplitude [mm]')
                    case 'imfdin'
                        ylabel('Amplitude [in]')
                    case 'imfdft'
                        ylabel('Amplitude [ft]')
                    case 'imfvm'
                        ylabel('Amplitude [m/s]')
                    case 'imfvmm'
                        ylabel('Amplitude [mm/s]')
                    case 'imfvin'
                        ylabel('Amplitude [in/s]')
                    case 'imfvft'
                        ylabel('Amplitude [ft/s]')
                    case 'imfam'
                        ylabel('Amplitude [m/s^2]')
                    case 'imfamm'
                        ylabel('Amplitude [mm/s^2]')
                    case 'imfain'
                        ylabel('Amplitude [in/s^2]')
                    case 'imfaft'
                        ylabel('Amplitude [ft/s^2]')
                    case 'imffn'
                        ylabel('Amplitude [N]')
                    case 'imffkn'
                        ylabel('Amplitude [kN]')
                    case 'imfflb'
                        ylabel('Amplitude [lb]')
                    case 'imffkip'
                        ylabel('Amplitude [kip]')
                    case 'phase'
                        ylabel('Phase [rad]')
                    case 'phasevel'
                        ylabel('Phase Velocity [rad/s]')
                    otherwise
                        if isempty(Label)
                            ylabel('Signal [\cdot]')
                        else
                            ylabel(Label)
                        end
                end
            else
                switch Label
                    case 'Acc'
                        ylabel({'Accelerance','[(m/s^2)/N]'})
                    case 'Vel'
                        ylabel({'Mobility','[(m/s)/N]'})
                    case 'Disp'
                        ylabel({'Receptance','[m/N]'})
                    case 'Dispmm'
                        ylabel({'Receptance','[mm/N]'})
                    case 'Tension'
                        ylabel('Tension/Input Force [N/N]')
                    otherwise
                        ylabel('')
                end
            end
        end


        function IMF = WBEMD(time,signal,minFreq,maxFreq,charFreq,options)
            % WBEMD: Extract monochromatic IMFs from provided time-series.
            %
            % Required Inputs
            % ---------------------------------------------
            % time - Time vector
            % signal - Column vector or matrix where each column represents
            %          one signal.
            % minFreq - Minimum frequency computed in the wavelet transform.
            % maxFreq - Maximum frequency computed in the wavelet transform.
            % charFreq – A cell that represents the characteristic frequencies
            %            of each signal that is being analyzed. Thus, if x is
            %            a p by 2 matrix (two signals are being analyzed)
            %            then Freq.FChar would be defined as
            %
            %            charFreq{1} = [Fn;Fn-1;...;F2;F1];
            %            charFreq{2} = [Gm;Gm-1;...;G2;G1];
            %
            %            where Fn,...,F1 and Gm,...,G1 are the characteristic
            %            time scales of the first and second time signals,
            %            respectively, and m and n are integers and can have
            %            different values.
            %
            % Optional Inputs (Uses Name-value format)
            % ---------------------------------------------
            % numFreq - The number of sample points contained in the interval
            %           [minFreq maxFreq]. The default value is 100, but
            %           better resolution in the frequency domain can be
            %           obtained by increasing this number.
            %
            % motherWaveletFreq - The frequency of the mother wavelet.
            %                     The default value is 2.
            %
            % MaskSig - Indicates whether the masking signal is constructed
            %           using a sine function or a cosine function. If the value
            %           is 1, then a sine is used. If the value is 2, then
            %           a cosine is used. The default value is 1.
            %
            % Outputs
            % ---------------------------------------------
            % IMF - A cell with the format: IMF{n}(:,m) which
            %       represents the mth IMF extracted from the nth signal
            %       decomposed by WBEMD. For example,
            %
            %       IMF{1}(:,1) is the 1st IMF extracted from the first
            %       signal decomposed by WBEMD.
            %
            %       IMF{4}(:,3) is the 3rd IMF extracted from for the
            %       4th signal decomposed by WBEMD.
            %
            % Example
            % ---------------------------------------------
            % clc,close all
            % Fs = 8192;
            % dt = 1/Fs;
            % time = 0:dt:10;
            % time = time';
            % signal = sin(10*2*pi*time)+sin(15*2*pi*time);
            % minFreq = 5;
            % maxFreq = 20;
            % charFreq{1} = [15;10];
            % Mirror.On = 0;
            % IMF = MoDAL.WBEMD(time,signal,minFreq,maxFreq,charFreq,"Mirror",Mirror);
            %
            % MoDAL.PlotTSWT_Compare(time,IMF{1}(:,2),time,IMF{1}(:,1),5,20,'legends',{'IMF 1','IMF 2'}
            %
            %
            % Dependencies
            % ---------------------------------------------
            % This code depends on the EMD package developed by Flandrin and Rilling
            % that can be downloaded from http://perso.ens-lyon.fr/patrick.flandrin/emd.html.
            % Most importantly, the code requires the mex version of emdc_fix that
            % comes with the above EMD package and it must be in the same folder as
            % wbemd. This requires that an appropriate compiler be installed on your
            % computer prior to the installation of the EMD package. A list of suitable
            % compilers is provided by Mathworks at https://www.mathworks.com/support/compilers.html.
            %
            % Copyright
            % ---------------------------------------------
            % Copyright (c) Keegan J. Moore, 2023
            % kmoore@unl.edu
            %
            % Citations
            % ---------------------------------------------
            % When using this code, please cite the following articles
            %
            % K.J. Moore, M. Kurt, M. Eriten, D.M. McFarland, L.A. Bergman, A.F.
            % Vakakis, "Wavelet-Bounded Empirical Mode Decomposition for Measured Time
            % Series Analysis," Mechanical Systems and Signal Processing, 99:14-29,
            % 2018. https://doi.org/10.1016/j.ymssp.2017.06.005
            %
            % K.J. Moore, M. Kurt, M. Eriten, D.M. McFarland, L.A. Bergman, A.F.
            % Vakakis, "Wavelet-Bounded Empirical Mode Decomposition for Vibro-Impact
            % Analysis," Nonlinear Dynamics, 93(3):1559-1577, 2018.
            % https://doi.org/10.1007/s11071-018-4276-0



            % Mirroring is not fully implemented and will be added in the
            % future. The list below is work in-progress.
            %
            % Mirror is a structure with the entries:
            %
            %   Mirror.On - Takes a value of 0 or 1, which indicate whether or not
            %               mirroring should be used in the decomposition.
            %
            %   Mirror.Freq - a cell that indicates which characteristic time
            %                scales mirroring will be applied to. For example,
            %
            %                Mirror.FMirror{1} = [Fn;F2;F1];
            %                Mirror.FMirror{2} = [Gm;Gm-1;G3;G1];
            %
            %                would result in mirroring being applied to the characteristic time
            %                scales Fn, F2, and F1 for the first signal and Gm, Gm-1, G3, and G1 for
            %                the second signal.
            %
            %   Mirror.Chop1 - A cell that indicates what percentage of the signal
            %                  should be mirrored at the beginning of the signal for a particular
            %                  characteristic frequency. The number of entries must equal the number
            %                  of entries in Mirror.FMirror and the recommended value is 0.2. For
            %                  example,
            %
            %                  Mirror.Chop1{1} = [0.2;0.2;0.2];
            %                  Mirror.Chop1{2} = [0.2;0.2;0.2;0.2];
            %
            %                  where the first and second lines correspond to the first and second
            %                  signals.
            %
            %   Mirror.Chop2 - The same as Mirror.Chop1, except that it applies to the
            %                  end of the signal instead of the beginning.
            %
            %
            % The default setting is for mirroring to be turned off and not applied to
            % the signal.

            arguments
                time double
                signal double
                minFreq double
                maxFreq double
                charFreq cell
                options.motherWaveletFreq double = 2;
                options.numFreq double = 100;
                options.maskSig double = 1;
                options.mirrorOn double = 0;
                options.mirrorFreq cell = {};
                options.mirrorChop1 cell = {};
                options.mirrorChop2 cell = {};
                options.Mirror struct = [];
            end
            if ~isfield(options.Mirror,'On'); options.Mirror(1).On = 0; end

            % Wavelet Parameters
            numFreq = options.numFreq;
            motherWaveletFreq = options.motherWaveletFreq;

            % Other Parameters
            breakers = 0;
            M = size(signal,2);

            EMD_OPTs.MAXMODES = 1;

            NumIMFs = 0;
            NumMaskSig = 0;
            for b = 1:M
                if b > 1
                    clear('X_emd','masksig')
                end

                if options.Mirror.On == 1
                    fmirr = options.Mirror.FMirror{b};
                    chops1 = options.Mirror.Chop1{b};
                    chops2 = options.Mirror.Chop2{b};
                end

                fg = charFreq{b};
                if size(time,2) == size(signal,2)
                    ti = time(:,b);
                else
                    ti = time;
                end
                X_emd = signal(:,b);
                L_FInt = length(fg);

                o = 1;
                m = 1;
                tic
                while(1)
                    if fg(o) < maxFreq/2.5
                        maxFreq = maxFreq/2;
                    end
                    if options.Mirror.On == 1
                        if fg(o) == fmirr(m)

                            chp1 = chops1(m);
                            chp2 = chops2(m);

                            L = length(ti); L_chp1 = MoDAL.chop(L*chp1,2); L_chp2 = MoDAL.chop(L*chp2,2);

                            [x_o_ini,~] = MoDAL.MirrorImgSigOIni(ti(1:L_chp1),X_emd(1:L_chp1));
                            [x_e_ini,~] = MoDAL.MirrorImgSigEIni(ti(1:L_chp1),X_emd(1:L_chp1));
                            [x_o_fin,~] = MoDAL.MirrorImgSigOFin(ti(end-L_chp2+1:end),X_emd(end-L_chp2+1:end));
                            [x_e_fin,~] = MoDAL.MirrorImgSigEFin(ti(end-L_chp2+1:end),X_emd(end-L_chp2+1:end));
                            if options.Mirror.Ini == 0
                                if options.Mirror.Fin == 0
                                    X_emd_mirror = [x_o_ini; X_emd; x_o_fin];

                                else
                                    X_emd_mirror = [x_o_ini; X_emd; x_e_fin];
                                end
                            else
                                if options.Mirror.Fin == 0
                                    X_emd_mirror = [x_e_ini; X_emd; x_o_fin];
                                else
                                    X_emd_mirror = [x_e_ini; X_emd; x_e_fin];
                                end
                            end
                            X_emd = X_emd_mirror;

                        end
                    end

                    if fg(o)/fg(o+1) > 20
                        fprintf('\nComputing IMF %g of of Time Series %g.\n',o,b)
                        F1 = fg(o);

                        if options.maskSig == 1
                            masksig = max(X_emd)*(sin(F1*ti));
                        elseif options.maskSig == 2
                            masksig = max(X_emd)*(cos(F1*ti));
                        end

                        if options.Mirror.On == 1
                            if fg(o) == fmirr(m)
                                chp1 = chops1(m);
                                chp2 = chops2(m);
                                L = length(ti); L_chp1 = MoDAL.chop(L*chp1,2); L_chp2 = MoDAL.chop(L*chp2,2);

                                [xm_o_ini, ~] = MoDAL.MirrorImgSigOIni(ti(1:L_chp1),masksig(1:L_chp1));
                                [xm_e_ini, ~] = MoDAL.MirrorImgSigEIni(ti(1:L_chp1),masksig(1:L_chp1));
                                [xm_o_fin, ~] = MoDAL.MirrorImgSigOFin(ti(end-L_chp2+1:end),masksig(end-L_chp2+1:end));
                                [xm_e_fin, ~] = MoDAL.MirrorImgSigEFin(ti(end-L_chp2+1:end),masksig(end-L_chp2+1:end));

                                if options.Mirror.Ini == 0
                                    if options.Mirror.Fin == 0
                                        masksig = [xm_o_ini;masksig;xm_o_fin];
                                    else
                                        masksig = [xm_o_ini;masksig;xm_e_fin];
                                    end
                                else
                                    if options.Mirror.Fin == 0
                                        masksig = [xm_e_ini;masksig;xm_o_fin];
                                    else
                                        masksig = [xm_e_ini;masksig;xm_e_fin];
                                    end
                                end
                            end
                        end

                        iter = 100;
                        [imf1,~] = emdc_fix(ti,X_emd+masksig,iter,1);
                        [imf2,~] = emdc_fix(ti,X_emd-masksig,iter,1);
                        if size(imf1,1) ~= size(imf2,1)
                            warning('emd:warning',['the two sets of IMFs have different sizes: ',int2str(size(imf1,1)),' and ',int2str(size(imf2,1)),' IMFs.'])
                        end
                        S1 = size(imf1,1);
                        S2 = size(imf2,1);
                        if S1 ~= S2
                            if S1 < S2
                                tmp = imf1;
                                imf1 = imf2;
                                imf2 = tmp;
                            end
                            imf2(max(S1,S2),1) = 0;
                        end
                        IMF_Proc(1,1,:,:) = (imf1+imf2)/2;
                        NumIMFs = NumIMFs+1;
                    else

                        NumMaskSig = NumMaskSig+1;
                        Opts = optimoptions('patternsearch','TolMesh',1e-8);
                        if size(fg,2) == 1
                            if o == 1
                                fr_std = abs(fg(o)-fg(o+1));
                            elseif abs(fg(o,1)-fg(o+1,1)) < abs(fg(o)-fg(o-1))
                                fr_std = abs(fg(o,1)-fg(o+1,1));
                            else
                                fr_std = abs(fg(o,1)-fg(o-1,1));
                            end
                        else
                            if (fg(o,1)-fg(o,2)) < (fg(o,2) - fg(o+1,1))
                                fr_std = fg(o,1)-fg(o,2);
                            else
                                fr_std = fg(o,2)-fg(o+1,1);
                            end
                        end

                        fprintf('Optimizing IMF %g of Time Series %g.\n ',o,b)

                        [freq,modx] = MoDAL.WaveletTransform(time,signal(:,b),minFreq,maxFreq, ...
                            'numFreq',numFreq,'motherWaveletFreq',motherWaveletFreq);
                        modx(isnan(modx))=0;
                        if options.Mirror.On == 1
                            if fg(o,1) == fmirr(m)
                                modx = modx(L_chp1:end-L_chp2+1,:);
                            end
                        end
                        maxmodx = max(modx);

                        freqx = sum(freq <= (fg(o,1) - (fg(o,1)-fg(o+1,1))/4));
                        modx_normer = max(maxmodx(freqx:end));
                        x0 = [1.1 1.1];
                        lb = [-1.3 1];
                        ub = [1.3 3];
                        if length(freq(maxmodx == modx_normer)) > 1
                            freq_normer = mean(freq(maxmodx == modx_normer));
                        else
                            freq_normer = freq(maxmodx == modx_normer);
                        end
                        check_hbound = 1.1;

                        ilen = 200;
                        jlen = 200;
                        mod_criteria1 = zeros(numFreq+1,ilen,jlen);
                        mod_criteria2 =zeros(numFreq+1,ilen,jlen);
                        checksval = zeros(ilen*jlen,1);
                        check_errtol = 0.01;
                        u = 1;
                        for j = 1:jlen
                            checkstd2 = 0.001;
                            checkstd1 = 3;
                            for i = 1:ilen
                                mod_criteria1(:,i,j) = normpdf(freq,freq_normer,fr_std/checkstd1)';
                                mod_criteria1(:,i,j) = mod_criteria1(:,i,j)/max(mod_criteria1(:,i,j))*check_hbound;
                                mod_criteria1(mod_criteria1(:,i,j) < check_errtol,i,j) = check_errtol;
                                checksval(u) = trapz(freq,mod_criteria1(:,i,j));
                                mod_criteria2(:,i,j) =  lognpdf(freq,log(freq_normer*1.0),checkstd2)';
                                mod_criteria2(:,i,j) = mod_criteria2(:,i,j)/max(mod_criteria2(:,i,j))*check_hbound;
                                mod_criteria2(mod_criteria2(:,i,j) < check_errtol,i,j) = check_errtol;
                                checkstd2 = checkstd2+0.001;
                                checkstd1 = checkstd1-0.01;
                                u = u+1;
                            end
                            check_errtol = check_errtol+0.001;
                        end

                        Omega = fg(o+1)+0.25*(fg(o)-fg(o+1));
                        gamma = 0.85;
                        delta = 0.25;


                        [xopt,~,~,~] = patternsearch(@(x)MoDAL.objfun(x,ti,X_emd,modx_normer,minFreq, ...
                            maxFreq,numFreq,motherWaveletFreq,mod_criteria1,mod_criteria2,checksval, ...
                            fg,gamma,delta,Omega),x0,[],[],[],[],lb,ub,[],Opts);

                        disp ' '
                        z = xopt;
                        if options.maskSig == 1
                            masksig = max(abs(X_emd))*(z(1)*sin(2*pi*fg(o)*z(2)*ti));
                        else
                            masksig = max(abs(X_emd))*(z(1)*cos(2*pi*fg(o)*z(2)*ti));
                        end

                        if options.Mirror.On == 1
                            if fg(o) == fmirr(m)
                                chp1 = chops1(m);
                                chp2 = chops2(m);
                                L = length(ti); L_chp1 = MoDAL.chop(L*chp1,2); L_chp2 = MoDAL.chop(L*chp2,2);

                                [xm_o_ini,~] = MoDAL.MirrorImgSigOIni(ti(1:L_chp1),masksig(1:L_chp1));
                                [xm_e_ini,~] = MoDAL.MirrorImgSigEIni(ti(1:L_chp1),masksig(1:L_chp1));
                                [xm_o_fin,~] = MoDAL.MirrorImgSigOFin(ti(end-L_chp2+1:end),masksig(end-L_chp2+1:end));
                                [xm_e_fin,~] = MoDAL.MirrorImgSigEFin(ti(end-L_chp2+1:end),masksig(end-L_chp2+1:end));

                                if options.Mirror.Ini == 0
                                    if Mirror.Fin == 0
                                        masksig = [xm_o_ini;masksig;xm_o_fin];
                                    else
                                        masksig = [xm_o_ini;masksig;xm_e_fin];
                                    end
                                else
                                    if options.Mirror.Fin == 0
                                        masksig = [xm_e_ini;masksig;xm_o_fin];
                                    else
                                        masksig = [xm_e_ini;masksig;xm_e_fin];
                                    end
                                end
                            end
                        end
                        try
                            [imf1,~] = emdc_fix(ti,X_emd+masksig,100,1);
                            [imf2,~] = emdc_fix(ti,X_emd-masksig,100,1);

                            S1 = size(imf1,1);
                            S2 = size(imf2,1);
                            if S1 ~= S2
                                if S1 < S2
                                    tmp = imf1;
                                    imf1 = imf2;
                                    imf2 = tmp;
                                end
                                imf2(max(S1,S2),1) = 0;
                            end
                            IMF_Proc = (imf1+imf2)/2;
                        catch
                            EMD_OPTs.MAXMODES = 1;
                            EMD_OPTs.MASK = masksig;
                            EMD_OPTs.MAXITERATIONS = 100;
                            [IMF_Proc,~] = emd(X_emd,EMD_OPTs);
                        end

                    end
                    toc
                    % This portion of code catches if an IMF cannot be extracted at
                    % a particular frequency and sets that IMF to zero and sets the
                    % remainder to the original time series. If the time signal is
                    % measured near a node of vibration the response at a
                    % particular frequency will be near zero and this algorithm
                    % will not be able to extract an IMF at that frequency. This
                    % code catches those occurances.
                    if breakers == 1
                        disp 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
                        fprintf('\nCheck failed on IMF %g with a frequency of %g of time series %g.\n',o,fg(o),b)
                        fprintf('Setting IMF to zero and proceeding.\n')
                        IMF_Proc(1,:) = zeros(size(squeeze(IMF_Proc(1,:))));
                        IMF_Proc(2,:) = X_emd';
                    end

                    if options.Mirror.On == 1
                        if fg(o,1) == fmirr(m)
                            IMF_Proc(1,:) = IMF_Proc(:,L_chp1:end-L_chp2+1);
                        end
                    end
                    IMF_Keep(o,:) = IMF_Proc(1,:);
                    Rem_Keep(o,:) = IMF_Proc(2,:);
                    X_emd = IMF_Proc(2,:)';
                    breakers = 0;
                    if options.Mirror.On == 1
                        if m == length(fmirr)
                            m = 1;
                        else
                            m = m+1;
                        end
                    end

                    if o >= L_FInt-1
                        break
                    end

                    o = o+1;
                    NumIMFs = NumIMFs+1;
                end
                EMD_OPTs.MASK = [];

                IMF{b} = [IMF_Keep;Rem_Keep(end,:)]';

                IMF_Keep = [];
                Rem_Keep = [];
            end
        end

        function [amp,theta,omega] = InstAmpFreq(time,signal,options)
            % Citation
            % --------
            % K.J. Moore, M. Kurt, M. Eriten, D.M. McFarland, L.A. Bergman, A.F. Vakakis,
            % “Wavelet-Bounded Empirical Mode Decomposition for Measured Time Series Analysis,”
            % Mechanical Systems and Signal Processing [IF: 8.934], 99:14–29, 2018.
            % https://dx.doi.org/10.1016/j.ymssp.2017.06.005
            arguments
                time (:,1) double
                signal (:,1) double
                options.cutoffFreq double = 1;
                options.filtOrder double = 3;
                options.mirrori string = 'e';
                options.mirrorf string = 'e';
                options.chp1 double = 0.2;
                options.chp2 double = 0.2;
            end
            dt = time(2)-time(1);

            [x_mirror,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd] = MoDAL.MirrorSignal(time,signal, ...
                options.mirrori,options.mirrorf,options.chp1,options.chp2);

            X = hilbert(x_mirror);

            theta1 = unwrap(angle(X));
            omega1 = gradient(theta1)/dt/(2*pi);
            amp1 = abs(X);
            Fs = 1/dt;
            Fnyq = Fs/2;

            Fcs = options.cutoffFreq/Fnyq;
            [b,a] = butter(options.filtOrder,Fcs,'low');
            theta2 = filtfilt(b,a,theta1);
            omega2 = filtfilt(b,a,omega1);
            amp2 = filtfilt(b,a,amp1);
            amp = MoDAL.UnmirrorSignal(amp2,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd);
            theta = MoDAL.UnmirrorSignal(theta2,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd);
            omega = MoDAL.UnmirrorSignal(omega2,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd);
        end

        function [fr,mods] = WaveletSignal(time,signal,minFreq,maxFreq,numFreq,motherWaveletFreq,mirrori,mirrorf)
            chp1 = 0.2;
            chp2 = 0.2;

            [x_mirror,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd] = MoDAL.MirrorSignal(time,signal,mirrori,mirrorf,chp1,chp2);

            [fr,mods] = MoDAL.WaveletTransform(time,x_mirror,minFreq,maxFreq, ...
                'motherWaveletFreq',motherWaveletFreq,'numFreq',numFreq);

            if NoMirrorIni == 0 && NoMirrorEnd == 0
                mods = mods(L_chp1:end-L_chp2+1,:);
            elseif NoMirrorIni == 0 && NoMirrorEnd == 1
                mods = mods(L_chp1:end);
            elseif NoMirrorIni == 1 && NoMirrorEnd == 0
                mods = mods(1:end-L_chp2+1);
            end
        end

    end

    properties (Constant,Access=private)
        UnbiasedOneMinusPinkMat = [   0.94107            1            1'
            0.93738      0.99608      0.99608;
            0.93369      0.99216      0.99216;
            0.93      0.98824      0.98824;
            0.92631      0.98431      0.98431;
            0.92262      0.98039      0.98039;
            0.91893      0.97647      0.97647;
            0.91524      0.97255      0.97255;
            0.91155      0.96863      0.96863;
            0.90786      0.96471      0.96471;
            0.90417      0.96078      0.96078;
            0.90047      0.95686      0.95686;
            0.89678      0.95294      0.95294;
            0.89309      0.94902      0.94902;
            0.8894       0.9451       0.9451;
            0.88571      0.94118      0.94118;
            0.88202      0.93725      0.93725;
            0.87833      0.93333      0.93333;
            0.87464      0.92941      0.92941;
            0.87095      0.92549      0.92549;
            0.86726      0.92157      0.92157;
            0.86357      0.91765      0.91765;
            0.85988      0.91373      0.91373;
            0.85619       0.9098       0.9098;
            0.8525      0.90588      0.90588;
            0.84881      0.90196      0.90196;
            0.84512      0.89804      0.89804;
            0.84143      0.89412      0.89412;
            0.83774       0.8902       0.8902;
            0.83405      0.88627      0.88627;
            0.83036      0.88235      0.88235;
            0.82667      0.87843      0.87843;
            0.82297      0.87451      0.87451;
            0.81928      0.87059      0.87059;
            0.81559      0.86667      0.86667;
            0.8119      0.86275      0.86275;
            0.80821      0.85882      0.85882;
            0.80452       0.8549       0.8549;
            0.80083      0.85098      0.85098;
            0.79714      0.84706      0.84706;
            0.79345      0.84314      0.84314;
            0.78976      0.83922      0.83922;
            0.78607      0.83529      0.83529;
            0.78238      0.83137      0.83137;
            0.77869      0.82745      0.82745;
            0.775      0.82353      0.82353;
            0.77131      0.81961      0.81961;
            0.76762      0.81569      0.81569;
            0.76393      0.81176      0.81176;
            0.76024      0.80784      0.80784;
            0.75655      0.80392      0.80392;
            0.75286          0.8          0.8;
            0.74917      0.79608      0.79608;
            0.74548      0.79216      0.79216;
            0.74178      0.78824      0.78824;
            0.73809      0.78431      0.78431;
            0.7344      0.78039      0.78039;
            0.73071      0.77647      0.77647;
            0.72702      0.77255      0.77255;
            0.72333      0.76863      0.76863;
            0.71964      0.76471      0.76471;
            0.71595      0.76078      0.76078;
            0.71226      0.75686      0.75686;
            0.70857      0.75294      0.75294;
            0.70488      0.74902      0.74902;
            0.70119       0.7451       0.7451;
            0.6975      0.74118      0.74118;
            0.69381      0.73725      0.73725;
            0.69012      0.73333      0.73333;
            0.68643      0.72941      0.72941;
            0.68274      0.72549      0.72549;
            0.67905      0.72157      0.72157;
            0.67536      0.71765      0.71765;
            0.67167      0.71373      0.71373;
            0.66798       0.7098       0.7098;
            0.66428      0.70588      0.70588;
            0.66059      0.70196      0.70196;
            0.6569      0.69804      0.69804;
            0.65321      0.69412      0.69412;
            0.64952       0.6902       0.6902;
            0.64583      0.68627      0.68627;
            0.64214      0.68235      0.68235;
            0.63845      0.67843      0.67843;
            0.63476      0.67451      0.67451;
            0.63107      0.67059      0.67059;
            0.62738      0.66667      0.66667;
            0.62369      0.66275      0.66275;
            0.62      0.65882      0.65882;
            0.61631       0.6549       0.6549;
            0.61262      0.65098      0.65098;
            0.60893      0.64706      0.64706;
            0.60524      0.64314      0.64314;
            0.60155      0.63922      0.63922;
            0.59786      0.63529      0.63529;
            0.59417      0.63137      0.63137;
            0.59048      0.62745      0.62745;
            0.58678      0.62353      0.62353;
            0.58309      0.61961      0.61961;
            0.5794      0.61569      0.61569;
            0.57571      0.61176      0.61176;
            0.57202      0.60784      0.60784;
            0.56833      0.60392      0.60392;
            0.56464          0.6          0.6;
            0.56095      0.59608      0.59608;
            0.55726      0.59216      0.59216;
            0.55357      0.58824      0.58824;
            0.54988      0.58431      0.58431;
            0.54619      0.58039      0.58039;
            0.5425      0.57647      0.57647;
            0.53881      0.57255      0.57255;
            0.53512      0.56863      0.56863;
            0.53143      0.56471      0.56471;
            0.52774      0.56078      0.56078;
            0.52405      0.55686      0.55686;
            0.52036      0.55294      0.55294;
            0.51667      0.54902      0.54902;
            0.51298       0.5451       0.5451;
            0.50928      0.54118      0.54118;
            0.50559      0.53725      0.53725;
            0.5019      0.53333      0.53333;
            0.49821      0.52941      0.52941;
            0.49452      0.52549      0.52549;
            0.49083      0.52157      0.52157;
            0.48714      0.51765      0.51765;
            0.48345      0.51373      0.51373;
            0.47976       0.5098       0.5098;
            0.47607      0.50588      0.50588;
            0.47238      0.50196      0.50196;
            0.46869      0.49804      0.49804;
            0.465      0.49412      0.49412;
            0.46131       0.4902       0.4902;
            0.45762      0.48627      0.48627;
            0.45393      0.48235      0.48235;
            0.45024      0.47843      0.47843;
            0.44655      0.47451      0.47451;
            0.44286      0.47059      0.47059;
            0.43917      0.46667      0.46667;
            0.43548      0.46275      0.46275;
            0.43179      0.45882      0.45882;
            0.42809       0.4549       0.4549;
            0.4244      0.45098      0.45098;
            0.42071      0.44706      0.44706;
            0.41702      0.44314      0.44314;
            0.41333      0.43922      0.43922;
            0.40964      0.43529      0.43529;
            0.40595      0.43137      0.43137;
            0.40226      0.42745      0.42745;
            0.39857      0.42353      0.42353;
            0.39488      0.41961      0.41961;
            0.39119      0.41569      0.41569;
            0.3875      0.41176      0.41176;
            0.38381      0.40784      0.40784;
            0.38012      0.40392      0.40392;
            0.37643          0.4          0.4;
            0.37274      0.39608      0.39608;
            0.36905      0.39216      0.39216;
            0.36536      0.38824      0.38824;
            0.36167      0.38431      0.38431;
            0.35798      0.38039      0.38039;
            0.35429      0.37647      0.37647;
            0.35059      0.37255      0.37255;
            0.3469      0.36863      0.36863;
            0.34321      0.36471      0.36471;
            0.33952      0.36078      0.36078;
            0.33583      0.35686      0.35686;
            0.33214      0.35294      0.35294;
            0.32845      0.34902      0.34902;
            0.32476       0.3451       0.3451;
            0.32107      0.34118      0.34118;
            0.31738      0.33725      0.33725;
            0.31369      0.33333      0.33333;
            0.31      0.32941      0.32941;
            0.30631      0.32549      0.32549;
            0.30262      0.32157      0.32157;
            0.29893      0.31765      0.31765;
            0.29524      0.31373      0.31373;
            0.29155       0.3098       0.3098;
            0.28786      0.30588      0.30588;
            0.28417      0.30196      0.30196;
            0.28048      0.29804      0.29804;
            0.27679      0.29412      0.29412;
            0.27309       0.2902       0.2902;
            0.2694      0.28627      0.28627;
            0.26571      0.28235      0.28235;
            0.26202      0.27843      0.27843;
            0.25833      0.27451      0.27451;
            0.25464      0.27059      0.27059;
            0.25095      0.26667      0.26667;
            0.24726      0.26275      0.26275;
            0.24357      0.25882      0.25882;
            0.23988       0.2549       0.2549;
            0.23619      0.25098      0.25098;
            0.2325      0.24706      0.24706;
            0.22881      0.24314      0.24314;
            0.22512      0.23922      0.23922;
            0.22143      0.23529      0.23529;
            0.21774      0.23137      0.23137;
            0.21405      0.22745      0.22745;
            0.21036      0.22353      0.22353;
            0.20667      0.21961      0.21961;
            0.20298      0.21569      0.21569;
            0.19929      0.21176      0.21176;
            0.19559      0.20784      0.20784;
            0.1919      0.20392      0.20392;
            0.18821          0.2          0.2;
            0.18452      0.19608      0.19608;
            0.18083      0.19216      0.19216;
            0.17714      0.18824      0.18824;
            0.17345      0.18431      0.18431;
            0.16976      0.18039      0.18039;
            0.16607      0.17647      0.17647;
            0.16238      0.17255      0.17255;
            0.15869      0.16863      0.16863;
            0.155      0.16471      0.16471;
            0.15131      0.16078      0.16078;
            0.14762      0.15686      0.15686;
            0.14393      0.15294      0.15294;
            0.14024      0.14902      0.14902;
            0.13655       0.1451       0.1451;
            0.13286      0.14118      0.14118;
            0.12917      0.13725      0.13725;
            0.12548      0.13333      0.13333;
            0.12179      0.12941      0.12941;
            0.1181      0.12549      0.12549;
            0.1144      0.12157      0.12157;
            0.11071      0.11765      0.11765;
            0.10702      0.11373      0.11373;
            0.10333       0.1098       0.1098;
            0.099643      0.10588      0.10588;
            0.095952      0.10196      0.10196;
            0.092262     0.098039     0.098039;
            0.088571     0.094118     0.094118;
            0.084881     0.090196     0.090196;
            0.08119     0.086275     0.086275;
            0.0775     0.082353     0.082353;
            0.073809     0.078431     0.078431;
            0.070119      0.07451      0.07451;
            0.066428     0.070588     0.070588;
            0.062738     0.066667     0.066667;
            0.059048     0.062745     0.062745;
            0.055357     0.058824     0.058824;
            0.051667     0.054902     0.054902;
            0.047976      0.05098      0.05098;
            0.044286     0.047059     0.047059;
            0.040595     0.043137     0.043137;
            0.036905     0.039216     0.039216;
            0.033214     0.035294     0.035294;
            0.029524     0.031373     0.031373;
            0.025833     0.027451     0.027451;
            0.022143     0.023529     0.023529;
            0.018452     0.019608     0.019608;
            0.014762     0.015686     0.015686;
            0.011071     0.011765     0.011765;
            0.0073809    0.0078431    0.0078431;
            0.0036905    0.0039216    0.0039216;
            0            0            0];
    end

    %%%%% Static Private methods
    methods (Static,Access=private)

        function TSPlot(time,signal,options)
            arguments
                time double
                signal double
                options.label string = '';
                options.color string = 'k';
                options.fontSize double = 12;
                options.timeStart double = time(1);
                options.timeEnd double = time(end);
                options.tsLim double = nan;
                options.linestyle string = '-';
            end
            plot(time,signal,options.color,'LineStyle',options.linestyle)
            MoDAL.TimeLabel(options.label)
            MoDAL.YLabel(options.label)
            set(gca,'FontSize',options.fontSize)
            xlim([options.timeStart options.timeEnd])
            if ~isnan(options.tsLim)
                if length(options.tsLim) == 1
                    ylim([-options.tsLim options.tsLim])
                else
                    ylim(options.tsLim)
                end
            end
        end

        function FTPlot(time,signal,minFreq,maxFreq,options)
            arguments
                time double
                signal double
                minFreq double
                maxFreq double
                options.label string = '';
                options.force double = [];
                options.color string = 'k';
                options.fontSize double = 12;
                options.linestyle string = '-';
            end
            dt = time(2)-time(1);
            freq = 0:1/time(end):1/(2*dt);
            L = length(time);
            if length(options.force) == L
                Fx = fft(signal);
                FF = fft(options.force);
                FFx = abs(Fx(1:length(freq))./FF(1:length(freq)));
                semilogy(freq,FFx,options.color,'LineStyle',options.linestyle)
                MoDAL.YLabel(options.label,1)
            else
                Fx = fft(signal);
                FFx = abs(Fx(1:length(freq)))*2/L;
                semilogy(freq,FFx,options.color,'LineStyle',options.linestyle)
                MoDAL.YLabel(options.label)
            end
            xlim([minFreq maxFreq])
            MoDAL.FreqLabel(options.label)
            set(gca,'FontSize',options.fontSize)
        end

        function MWTPlot(freq,maxmods,minFreq,maxFreq,options)
            arguments
                freq double
                maxmods double
                minFreq double
                maxFreq double
                options.label string = '';
                options.color string = 'k';
                options.fontSize double = 12;
            end
            plot(freq,maxmods,options.color)
            xlim([minFreq maxFreq])
            MoDAL.FreqLabel(options.label)
            ylabel('Norm. Ampl. [\cdot]')
            set(gca,'FontSize',options.fontSize)
        end

        function WTSpectraPlot(time,freq,mods,options)
            imagesc(time,freq,(mods').^options.wtPower)
            set(gca,'ydir','nor')
            colormap(options.colorMap)
            MoDAL.WTLabel(options.label)
            xlim([options.timeStart options.timeEnd])
            set(gca,'FontSize',options.fontSize)
        end

        function TimeLabel(label)
            if contains(label,'ND')
                xlabel('Time [\cdot]')
            else
                xlabel('Time [s]')
            end
        end

        function WTLabel(label)
            if contains(label,'ND')
                xlabel('Time [\cdot]')
                ylabel('Freq. [\cdot]')
            else
                xlabel('Time [s]')
                ylabel('Freq. [Hz]')
            end
        end

        function FreqLabel(label)
            if contains(label,'ND')
                xlabel('Frequency [\cdot]')
            else
                xlabel('Frequency [Hz]')
            end
        end

        function P = objfun(z,time,X_emd,modx_normer,minFreq,maxFreq,numFreq,motherWaveletFreq,mod_criteria1,mod_criteria2,checksval,fg,gamma,delta,Omega)
            masksig = max(abs(X_emd))*(z(1)*sin(2*pi*fg(1)*z(2)*time));
            [imf1,~] = emdc_fix(time,X_emd+masksig,100,1);
            [imf2,~] = emdc_fix(time,X_emd-masksig,100,1);
            S1 = size(imf1,1);
            S2 = size(imf2,1);
            if S1 ~= S2
                if S1 < S2
                    tmp = imf1;
                    imf1 = imf2;
                    imf2 = tmp;
                end
                imf2(max(S1,S2),1) = 0;
            end
            IMF = (imf1+imf2)/2;

            [fr,modIMF] = MoDAL.WaveletTransform(time,IMF(1,:)',minFreq,maxFreq, ...
                'numFreq',numFreq,'motherWaveletFreq',motherWaveletFreq);

            modIMF(isnan(modIMF)) = 0;
            maxmodIMF = max(modIMF)/modx_normer;

            if max(maxmodIMF) < gamma || max(maxmodIMF) > 1.25
                P = 1e15;
            elseif max(maxmodIMF(fr <= fg(2))) > delta
                P = 1e15;
            elseif fr(maxmodIMF == max(maxmodIMF)) < Omega
                P = 1e15;
            else
                checkstest = bsxfun(@minus,mod_criteria1,maxmodIMF');
                checkstest2 = bsxfun(@minus,mod_criteria2,maxmodIMF');
                chs(:,:) = sum(checkstest < 0);
                chs2(:,:) = sum(checkstest2 < 0);
                [minval,minIndexes1] = min(chs(:));
                [minval2,~] = min(chs2(:));

                if minval ~= 0
                    if minval2 ~= 0
                        P = checksval(end)*2;
                    else
                        P = checksval(end);
                    end
                elseif minval2 ~= 0
                    P = checksval(end);
                else
                    P = checksval(minIndexes1);
                end
            end
        end

        function [x_mirror,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd] = ...
                MirrorSignal(time,signal,mirrori,mirrorf,chp1,chp2)
            arguments
                time (:,1) double
                signal (:,1) double
                mirrori string = 'e';
                mirrorf string  = 'e';
                chp1 double = 0.2;
                chp2 double = 0.2;
            end
            L = length(time); L_chp1 = MoDAL.chop(L*chp1,2); L_chp2 = MoDAL.chop(L*chp2,2);
            [x_o_ini,~] = MoDAL.MirrorImgSigOIni(time(1:L_chp1),signal(1:L_chp1));
            [x_e_ini,~] = MoDAL.MirrorImgSigEIni(time(1:L_chp1),signal(1:L_chp1));
            [x_o_fin,~] = MoDAL.MirrorImgSigOFin(time(end-L_chp2+1:end),signal(end-L_chp2+1:end));
            [x_e_fin,~] = MoDAL.MirrorImgSigEFin(time(end-L_chp2+1:end),signal(end-L_chp2+1:end));

            NoMirrorIni = 0;
            NoMirrorEnd = 0;
            if strcmp(mirrori,'e')
                if strcmp(mirrorf,'e')
                    x_mirror = [x_e_ini; signal; x_e_fin];
                elseif strcmp(mirrorf,'o')
                    x_mirror = [x_e_ini; signal; x_o_fin];
                else
                    x_mirror = [x_e_ini; signal];
                    NoMirrorEnd = 1;
                end
            elseif strcmp(mirrori,'o')
                if strcmp(mirrorf,'e')
                    x_mirror = [x_o_ini; signal; x_e_fin];
                elseif strcmp(mirrorf,'o')
                    x_mirror = [x_o_ini; signal; x_o_fin];
                else
                    x_mirror = [x_o_ini; signal];
                    NoMirrorEnd = 1;
                end
            else
                if strcmp(mirrorf,'e')
                    x_mirror = [signal; x_e_fin];
                elseif strcmp(mirrorf,'o')
                    x_mirror = [signal; x_o_fin];
                else
                    x_mirror = signal;
                    NoMirrorEnd = 1;
                end
                NoMirrorIni = 1;
            end
        end

        function Signal = UnmirrorSignal(signal,L_chp1,L_chp2,NoMirrorIni,NoMirrorEnd)
            if NoMirrorIni == 0 && NoMirrorEnd == 0
                Signal = signal(L_chp1:end-L_chp2+1);
            elseif NoMirrorIni == 0 && NoMirrorEnd == 1
                Signal = signal(L_chp1:end);
            elseif NoMirrorIni == 1 && NoMirrorEnd == 0
                Signal = signal(1:end-L_chp2+1);
            else
                Signal = signal;
            end
        end

        function X = chop(Xin,n,unit)
            %CHOP   CHOP(X,n) rounds elements of X to n significant figures.
            %       CHOP(X,n,unit) rounds the elements of X to n significant
            %   figures whose digits (mantissa) are exactly divisible
            %   by unit.
            %
            %       e.g. chop(3.141592,5)   returns 3.141600000..
            %       e.g. chop(3.141592,3,5) returns 3.150000000..
            %            chop(3.141592,3,3) returns 3.150000000..
            %            chop(3.141592,3,2) returns 3.140000000..
            %

            %   Copyright 1986-2002 The MathWorks, Inc.

            % Set last sig. fig. rounding to 1 if only two input arguments.
            if nargin<3
                unit=1;
            end

            % Cater for -ve numbers  and numbers = 0.
            X = abs(Xin) +(Xin==0);
            [nx,mx] = size(X);
            exponent = unit.*((10*ones(nx,mx)).^(floor(log10(X))-n+1));
            X = round(X./exponent).*exponent;

            % Put back sign and zeros
            X = sign(Xin).*X.*(Xin~=0);

        end

        function y = normpdf(x,mu,sigma)
            %NORMPDF Normal probability density function (pdf).
            %   Y = NORMPDF(X,MU,SIGMA) returns the pdf of the normal distribution with
            %   mean MU and standard deviation SIGMA, evaluated at the values in X.
            %   The size of Y is the common size of the input arguments.  A scalar
            %   input functions as a constant matrix of the same size as the other
            %   inputs.
            %
            %   Default values for MU and SIGMA are 0 and 1 respectively.
            %
            %   See also NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMRND, NORMSTAT.
            %
            %   References:
            %      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
            %          Distributions, 2nd ed., Wiley, 170pp.
            %
            %   Copyright 1993-2004 The MathWorks, Inc.


            if nargin<1
                error(message('stats:normpdf:TooFewInputs'));
            end
            if nargin < 2
                mu = 0;
            end
            if nargin < 3
                sigma = 1;
            end

            % Return NaN for out of range parameters.
            sigma(sigma <= 0) = NaN;

            try
                y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
            catch
                error(message('stats:normpdf:InputSizeMismatch'));
            end
        end

        function HideX
                % Extract axes from current figure
                Ax = gcf().Children(strcmp(get(gcf().Children,'type'),'axes'));
                switch length(Ax)
                    case 3
                        Ax3 = Ax(3);
                        Ax2 = Ax(2);
                        Ax1 = Ax(1);
                        Ax3.XTickLabel = {''};
                        Ax3.XLabel = xlabel('');
                        Ax3.Position(4) = 0.24;
                        Ax2.XTickLabel = {''};
                        Ax2.XLabel = xlabel('');
                        Ax2.Position(4) = 0.24;
                        Ax1.Position(4) = 0.24;
                        Ax1.XLabel = xlabel('Time [s]');
                    case 2
                        Ax(2).XTickLabel = '';
                        Ax(2).XLabel = xlabel('');
                end
        end

        function [x_mirror, t_mirror] = MirrorImgSigOIni(t, x)
            % Generate the mirror image signal by an odd symmetry about t=t0
            L = length(t); xtmp = x - x(1);
            t_mirror = [-(t(L:-1:2)-t(1))+t(1)];
            x_mirror = [-xtmp(L:-1:2)+x(1)];
        end

        function [x_mirror, t_mirror] = MirrorImgSigEIni(t, x)
            % Generate the mirror image signal by an even symmetry about t=t0
            L = length(t); xtmp = x - x(1);
            t_mirror = [-(t(L:-1:2)-t(1))+t(1)];
            x_mirror = [xtmp(L:-1:2)+x(1)];
        end

        function [x_mirror, t_mirror] = MirrorImgSigOFin(t, x)
            % Generate the mirror image signal by an odd symmetry about t=tf
            dt = t(2)-t(1); L = length(t); xtmp = x - x(end);
            t_mirror = [t(end)+dt*[1:L-1]'];
            x_mirror = [-xtmp(end-1:-1:end-L+1)+x(end)];
        end

        function [x_mirror, t_mirror] = MirrorImgSigEFin(t, x)
            % Generate the mirror image signal by an even symmetry about t=tf
            dt = t(2)-t(1); L = length(t); xtmp = x - x(end);
            t_mirror = [t(end)+dt*[1:L-1]'];
            x_mirror = [xtmp(end-1:-1:end-L+1)+x(end)];
        end
    end
end