classdef saObj < handle
    %UNTITLED2 Summary of this class goes here

    % setaddress (address)
    % ini (void): connect
    % config1 (void): configuration 1 for RE
    % config2 (void): configuration 1 for CE
    % config3 (void): configuration 1 for OATS?
    %   -import(id) :
    %                 id='limitRE'|'limitCE'|'limitOATs'|'corrRE'|'corrCE'|'corrOATs'|
    %                 load all calibration date from the MXE local drive
    %                 c:\data
    % meas   (void) : scan, search, create signal list and measure. 
    % remeas (void) : re-measure the list
    % takedata (void):  restart to make measurement.  
    % plot(void):     take screen shot 
    % 
    % scan
    % signal
    properties
        %address='USB0::0x0957::0x0f0b::MY52130013::0::INSTR';
        %address='USB0::0x0957::0x0F0B::MY51210168::0::INSTR';
        address='USB0::0x2A8D::0x0F0B::MY56400038::0::INSTR';
        sa=[];
        status='Closed';
        data=[]; % store the table scan value one trace only. 
        dataH=[]; % store the H scan value one trace only. 
        count=[];
        freq=[];
        freqH=[];
        RawList=[];
        SelectedList=[];
        peaktable;
        peakresult;
        opc_status;
        tw_location_for_peak; % single value. holder for the antenna postiion for the maximumum peak. 
        
    end
    
    methods
          function s=saObj
          end
        function setaddress(s,a)
            s.address=a;
        end
        function fs=fstart(s)
            fprintf(s.sa,'FREQ:STAR?')
            fs=str2num(fscanf(s.sa));
        end
        function fe=fstop(s)
            fprintf(s.sa,'FREQ:STOP?')
            fe=str2num(fscanf(s.sa));
        end
        
        function ini(s)
            s.sa = visa('agilent',s.address);
            % ni means National Instruments
            % open the session object 
            pause(0.1);
            s.sa.InputBufferSize=512000;
            fopen(s.sa);
            s.status= get(s.sa, 'Status')
        end
        function configCE(s)
            % Test if it is open
            fprintf(s.sa,'*RST;*CLS;*WAI\n');
            fprintf(s.sa,'INST:SEL EMI');
            fprintf(s.sa,':SYSTem:PRESet');
              while ~s.opc
            end
            %returns the measreument local variables in the current
            %measruement to their preset values
            fprintf(s.sa,'CONF:FSC');
            %s.import('scanTableCE');
            s.import('limitCE');
            s.import('corrCE');
            % Turn off built in amplifier, Noise floor remain the same
            fprintf(s.sa,':FSCan:SCAN5:POWer:GAIN 0');
            % set RF input port to 2
            fprintf(s.sa,':FSC:SCAN2:FEED:RF:PORT RFIN2');
            fprintf(s.sa,':FSCan:SCAN5:STATe 0');
            fprintf(s.sa,':FSCan:SCAN2:STATe 1');

            % set the sweep to continious
            fprintf(s.sa,'INITiate2:CONTinuous 1');
             %enable meter display
            fprintf(s.sa,'DISP:MET1 1');
            fprintf(s.sa,'DISP:MET2 1');
            fprintf(s.sa,'DISP:MET3 1');
           %set to single mode first
           fprintf(s.sa,'CALC:FSC:SLIS:DEL:ALL');
           %set to running mode
           fprintf(s.sa,'INIT2:RESTart');
            
        end
         function configRE(s)
            fprintf(s.sa,'*RST;*CLS;*WAI\n');
            fprintf(s.sa,'INST:SEL EMI');
            while ~s.opc
            end
            fprintf(s.sa,':SYSTem:PRESet');
            while ~s.opc
            end
          
  %          fprintf(s.sa,'CONF:FSC');
           % fprintf(s.sa,'INITiate:FSC');
           % s.import('scanTableRE');
            s.import('limitRE');
            s.import('corrRE');
            % Turn on built in amplifier, Noise floor drop 10dB
            fprintf(s.sa,':FSCan:SCAN5:POWer:GAIN 1');
              % set RF input port to 1
            fprintf(s.sa,':FSC:SCAN2:FEED:RF:PORT RFIN');
           % select scantable range
            fprintf(s.sa,':FSCan:SCAN5:STAR 30e6'); 
            fprintf(s.sa,':FSCan:SCAN5:STATe 1');
            fprintf(s.sa,':FSCan:SCAN2:STATe 0');
          
           % set the sweep to continious
            fprintf(s.sa,'INITiate2:CONTinuous 1');
            %enable meter display
            fprintf(s.sa,'DISP:MET1 1');
            fprintf(s.sa,'DISP:MET2 1');
            fprintf(s.sa,'DISP:MET3 1');
           %set to single mode first
           fprintf(s.sa,'CALC:FSC:SLIS:DEL:ALL');
            %set to running mode
           fprintf(s.sa,'INIT2:RESTart');
           
            
         end
         function configOATs(s)
            fprintf(s.sa,'*RST;*CLS;*WAI\n');
            fprintf(s.sa,'INST:SEL EMI');
            fprintf(s.sa,':SYSTem:PRESet');
         
          %  s.import('scanTableOATs');
            s.import('limitOATs');
            s.import('corrOATs');
           % Turn on built in amplifier, Noise floor drop 10dB
            fprintf(s.sa,':FSCan:SCAN5:POWer:GAIN 1');
            % set RF input port to 1
            fprintf(s.sa,':FSC:SCAN2:FEED:RF:PORT RFIN');
            % select scantable range;
            fprintf(s.sa,':FSCan:SCAN5:STATe 1');
            fprintf(s.sa,':FSCan:SCAN2:STATe 0');
 
            % set the sweep to continious
            fprintf(s.sa,'INITiate2:CONTinuous 1');
           %set to single mode first
           fprintf(s.sa,'CALC:FSC:SLIS:DEL:ALL');
            %set to running mode
           fprintf(s.sa,'INIT2:RESTart');
            
         end
          function StripChartMode(s)
           fprintf(s.sa,'DISP:SCH:VIEW:WIND:TRAC:X:MAX:DUR 10s');
           fprintf(s.sa,'CONFigure:SCHart');
       
           % fprintf(s.sa,'DISP:SCH:VIEW:WIND:TRAC:X:MAX:FULL');
            fprintf(s.sa,'INITiate:SCHart');
            fprintf(s.sa,'SCH:ABOR');
            
            fprintf(s.sa,'INIT:REST');
            fprintf(s.sa,'FORMat:DATA REAL,32\n') ;
            fprintf(s.sa,'FORMat:BORDer SWAPped');
            fprintf(s.sa,'FETCH:SCHart3?');
           [data, count] = binblockread(s.sa, 'float32');
           
            % set RF input port to 1
            
            
          end
           function opc_status=opc(s)
            fprintf(s.sa,'*OPC?');
            opc_status=str2num(strtrim(fscanf(s.sa)));
            s.opc_status=opc_status;
           end
            
          function Config_remeas_start(s,i,type) 
          % set span to 5MHz half the RBW
             % set the meter to SLIST number mm
          switch type
              case {'RE','OAT'}
                                  fprintf(s.sa,['CALC:FSC:SLIS:SET:MET ' num2str(i)]); % remeasure curretnt selected signal           
                              while ~s.opc
                             end
                            % move maker to the meter freq
                            fprintf(s.sa,['CALC:FSC:MARK1:MET ' num2str(i)]); % r       
                             while ~s.opc
                            end
                          s.setSAMode;   

                          while ~s.opc
                          end
                     %      fprintf(s.sa,'INIT2:CONT ON'); % set trigger mode into continius
                            fprintf(s.sa,'DISP:WIND:TRAC:Y:RLEV 80DBUVM'); % set trigger mode into single 
                            fprintf(s.sa,'FORMat:DATA INT,32\n') ;
                            fprintf(s.sa,'FORMat:BORDer SWAPped');
                            fprintf(s.sa,'INIT:CONT OFF'); % set trigger mode into continius
              case {'CE'}
                       fprintf(s.sa,['CALC:FSC:SLIS:SET:MET ' num2str(i)]); % remeasure curretnt selected signal           
                              while ~s.opc
                             end
                            % move maker to the meter freq
                            fprintf(s.sa,['CALC:FSC:MARK1:MET ' num2str(i)]); % r       
                             while ~s.opc
                            end
          end
          
                         
          end
        
           function Config_remeas_stop(s) 
          % set span to 5MHz half the RBW
         
              fprintf(s.sa,'INST:SEL EMI');
              while ~s.opc
              end
              
              fprintf(s.sa,'INIT:CONT OFF'); % set trigger mode into single 
                  end
                  
                  
           function data1=takedata(s,n) % which trace to take
           if ~strcmp(s.status,'open')
                msgbox('Spectrum Analyzer is not ready!')
           else
               fprintf(s.sa,'FORMat:DATA INT,32\n') ;
               fprintf(s.sa,'FORMat:BORDer SWAPped');
               fprintf(s.sa,['TRAC:FSC? TRACE' num2str(n)]);
               %fprintf(sa,'TRAC:DATA TRACE1');
               [data, count] = binblockread(s.sa, 'int32');
               data2=data*1e-3+107;
               s.count=count;
               data1=floor(data2.*10)./10;
               s.data=data1;
               
           end
           
        end
        function scanH(s,n) % which trace to take
           if ~strcmp(s.status,'open')
                msgbox('Spectrum Analyzer is not ready!')
           else
               fprintf(s.sa,'INIT:IMM'); % set trigger mode into single
              while ~s.opc
              end
              
             %  fprintf(s.sa,['TRAC:FSC? TRACE' num2str(n)]);
              % fprintf(s.sa,['READ:SAN' num2str(n) '?']);
               fprintf(s.sa,'TRAC? TRACE1');
              [data, count] = binblockread(s.sa, 'int32');
               while ~s.opc
               end              
               s.dataH=data*1e-3+107;
           end
               s.dataH=floor(s.dataH.*10)./10;
           
        end
        
        function takefreq(s)
               fprintf(s.sa,':FREQ:STAR?');
               fstart=str2num(fscanf(s.sa));
               fprintf(s.sa,':FREQ:STOP?');
               fstop=str2num(fscanf(s.sa));
               mm=size(s.data,1);
               df=(fstop-fstart)/(mm-1);
               freq=[fstart:df:fstop];
               s.freq=freq;
               
              
               
        end
        
        function search(s) % scan, search and no measure signallist. 
           %set to single mode first
           fprintf(s.sa,'CALC:FSC:SLIS:DEL:ALL');
           fprintf(s.sa,'INITiate2:CONTinuous 0');
           %scan
           fprintf(s.sa,':FSC:SEQ SCAN');
           fprintf(s.sa,'INIT2:IMM;*WAI;');
           fprintf(s.sa,'CALC:FSC:SLIS:DEL:ALL');
            %enable meter display
            fprintf(s.sa,'DISP:MET1 1');
            fprintf(s.sa,'DISP:MET2 1');
            fprintf(s.sa,'DISP:MET3 1');
           % search and measure
                        
            %fprintf(s.sa,':FSC:SEQ SAMeasure'); % search measure
           fprintf(s.sa,':FSC:SEQ SEAR');
           fprintf(s.sa,'INIT2:IMM;*WAI;'); % execute at read step.
       %    fprintf(s.sa,'FORMat:TRACe INT,32\n') ;
           s.readSLIST;
        end
        function Config_remeasure(s)
          % inittilize the MXXE
            s.configRE;
            fprintf(s.sa,'INIT2:CONT OFF'); % set trigger mode into single 
            fprintf(s.sa,':FSC:SEQ:REMeasure CURR'); % remeasure curretnt selected signal 
               % set global center frquency from EMI mode to SA mode on 
            fprintf(s.sa,'GLOBal:FREQuency:CENTer ON');
            while ~s.opc
            end
            % write the peaktable into the SLIS
            s.writeSLIST
        end
        
        function Config_remeasureCE(s)
            % set the MXE to remeasure CE signal list. 
          % inittilize the MXXE
            s.configCE;
            fprintf(s.sa,'INIT2:CONT OFF'); % set trigger mode into single 
            fprintf(s.sa,':FSC:SEQ:REMeasure CURR'); % remeasure curretnt selected signal 
               % set global center frquency from EMI mode to SA mode on 
            fprintf(s.sa,'GLOBal:FREQuency:CENTer ON');
            while ~s.opc
            end
            % write the peaktable into the SLIS
            s.writeSLIST
        end
        
        function remeasure(s)
            
            fprintf(s.sa,'INIT2:CONT OFF'); % set trigger mode into single 
            fprintf(s.sa,'CONFigure:EMI:SLISt CURR');
            fprintf(s.sa,':FSC:SEQ:REMeasure CURR'); % remeasure curretnt selected signal 
          %  fprintf(s.sa,['CALC:FSC:SLIS:SET:MET ' num2str(id)]); % remeasure curretnt selected signal    
         %  fprintf(s.sa,['CALC:FSC:SLIS:SNAP:MET']);
         %   fprintf(s.sa,'INIT2:IMM '); % remeasure curretnt selected signal           
            fprintf(s.sa,'INIT2:REST;*WAI'); % remeasure curretnt selected signal           
            
%                % set particular freq to remeasure
%                fprintf(s.sa,['CALC:FSC:SLIS:MARK:SIGN ' num2str(id)]);
%                 fprintf(s.sa,':FSC:SEQ:REMeasure MARK'); % remeasure at mark only 
          end
          function setsweep2(s,y)
          if y
                 fprintf(s.sa,'INIT2:CONT ON'); % set trigger mode into single 
          else
                 fprintf(s.sa,'INIT2:CONT OFF'); % set trigger mode into single 
          end
          
          end
          
          function setsweep(s,y)
          if y
                 fprintf(s.sa,'INIT:CONT ON'); % set trigger mode into single 
          else
                 fprintf(s.sa,'INIT:CONT OFF'); % set trigger mode into single 
          end
          
          end
          function setSAMode(s)
          
          fprintf(s.sa,'INST:SEL SA');
          while ~s.opc
          end
          fprintf(s.sa,'FREQ:SPAN 5000000');
          fprintf(s.sa,'UNIT:POW DBUVM'); % set the unit to dBuV
          
          end
          function setcenterfreq(s,y)
              
          
          fprintf(s.sa,['FREQ:CENTer ' num2str(y)]); % set the unit to dBuV
          
          end
          
        function readSLIST(s)
           fprintf(s.sa,'FORMat:DATA REAL,32\n') ;
           % fprintf(s.sa,'FORMat:BORDer NORMal');
           fprintf(s.sa,'FORMat:BORDer SWAPped');
           fprintf(s.sa,'FETCH:FSCan1?');
           [data, count] = binblockread(s.sa, 'float32');
            while ~s.opc
            end
           if count<9
               s.peakresult=[];
            else
               k=floor((length(data)-1)/9);
               for i=0:k-1
                   s.peakresult(i+1,:)=data(i*9+2:i*9+10);
               end
           end 
            s.peakresult=floor(s.peakresult.*10)./10; % format
        end
        
        function writeSLIST(s)
      %
       %     fprintf(s.sa,'CALC:FSC:MARK:MODE POS');
       % delete all exisiting signal list
             fprintf(s.sa,'CALC:FSC:SLIS:DEL:ALL');
             %enable marker
              fprintf(s.sa,'CALC:FSC:MARK:MODE1 POSition');
             mm=size(s.peaktable,1);
             for i=1:mm
              fprintf(s.sa,['CALC:FSC:MARK1:X ' num2str(s.peaktable.freq(i))]);
              while ~s.opc
               end
              fprintf(s.sa,'CALC:FSC:MARK1:SLIS');
                  while ~s.opc
               end
             end
             
        end
        function replace_peaktable(s,i)
            % delete the current signal on the signal list
            %fprintf(s.sa,['CALC:FSC:SLIS:DEL:SIGN ' num2str(i)]);
            %move the marker to the current new selected freq
            fprintf(s.sa,['CALC:FSC:MARK1:X ' num2str(s.peaktable.freq(i))]);
            %move the detector to the marker
            %fprintf(s.sa,['CALC:FSC:MARK1:TO:MET']);
            %snap signal to the cloest meter
            fprintf(s.sa,['FREQ:CENT ' num2str(s.peaktable.freq(i))]);
            fprintf(s.sa,['CALC:FSC:SLIS:REPL:MET ' num2str(i)]);
        %    fprintf(s.sa,['CALC:FSC:SLIS:SNAP:MET']);
            
                while ~s.opc
                end
%             %move the current maker freq into the signal list
%              fprintf(s.sa,'CALC:FSC:MARK1:SLIS');
%             % sort the freq. 
%             fprintf(s.sa,'CALC:FSC:SLIS:SORT:TYPE FREQ');
           
           
              
        end
        
        
        function selectSLIST(s)
            s.SelectedList=s.RawList;
            s.writeSLIST;
        end
        
        
        
        
         function plot(s)
             figure
         semilogx(s.freq,s.data) 
         xlabel('Frequency Hz');ylabel('dBuV');
         grid on
         
         
         end 
         function import(s,id)
         switch id
             case 'limitRE'
                 % HAS TO BE SAVED TO MYDOCMENT\SA\DATA\LIMITS
                  fprintf(s.sa,':MMEMory:LOAD:LIM LLINE1,"d:\data\RE3M.lim"');
                  fprintf(s.sa,':MMEMory:LOAD:LIM LLINE2,"d:\data\RE3MM.lim"');
                  fprintf(s.sa,'CALC:FSC:LLIN1:DISP ON');
                  fprintf(s.sa,'CALC:FSC:LLIN2:DISP ON');
                   
             case 'limitCE'
                  fprintf(s.sa,':MMEMory:LOAD:LIM LLINE3,"d:\data\CEQ.lim"');   
                  fprintf(s.sa,':MMEMory:LOAD:LIM LLINE4,"d:\data\CEA.lim"');
                  fprintf(s.sa,'CALC:FSC:LLIN3:DISP ON');
                  fprintf(s.sa,'CALC:FSC:LLIN4:DISP ON');
             case 'limitOATs'
                  fprintf(s.sa,':MMEMory:LOAD:LIM LLINE1,"d:\data\REOATS.lim"') ; 
                  fprintf(s.sa,':MMEMory:LOAD:LIM LLINE2,"d:\data\REOATSM.lim"'); 
                  fprintf(s.sa,'CALC:FSC:LLIN5:DISP ON');
                  fprintf(s.sa,'CALC:FSC:LLIN6:DISP ON');
                  
             case 'trace'
                  fprintf(s.sa,'MMEMory:LOAD:TRACE:DATA TRACE1 "myTrace.lim"');
             case 'corrRE'
                  fprintf(s.sa,'MMEMory:LOAD:CORRection 1,"c:\data\corrRE.csv"');
                  fprintf(s.sa,'SENS:CORR:CSET1 1');
                  fprintf(s.sa,':DISPlay:FSCan:VIEW:WIND:TRAC:Y:RLEVel 87dBuV');
                  
                  
                  % file should be pre-saved at c:\data\
             case 'corrCE'
                  fprintf(s.sa,'MMEMory:LOAD:CORRection 1,"c:\data\corrCE.csv"');
                  fprintf(s.sa,'SENS:CORR:CSET1 1');
                  fprintf(s.sa,':DISPlay:FSCan:VIEW:WIND:TRAC:Y:RLEVel 87dBuV');                  
                  % file should be pre-saved at c:\data\
             case 'corrOATs'
                  fprintf(s.sa,'MMEMory:LOAD:CORRection 1,"c:\data\corrOATs.csv"');
                  fprintf(s.sa,'SENS:CORR:CSET1 1');
                  fprintf(s.sa,':DISPlay:FSCan:VIEW:WIND:TRAC:Y:RLEVel 87dBuV');
                  % file should be pre-saved at c:\data\ 
             case 'signallist'
                  fprintf(s.sa,'MMEMory:LOAD:SLIST "mySignalList.csv"');
             case 'scantable'
                  fprintf(s.sa,'MMEMory:LOAD:SCAN "myScanTable.csv"');
             case 'setup'
                  fprintf(s.sa,':MMEMory:LOAD:STAT 1,"C:\DTDFILES\STATE000.STA"');
                  %*WAI\n
                  
         end
         
         end
         
    end
    
end


% %
%             fprintf(s.sa,'DISP:WIND:TRAC:Y:RLEV 86 DBUV')
%             fprintf(s.sa,'UNIT:POW DBUV')
%             fprintf(s.sa,'DISP:WIND:TRAC:Y:PDIV 10 DB')
% 
%             % [:SENSe]:DETector[:FUNCtion] NORMal|AVERage|POSitive|SAMPle|NEGative|QPEak
%             % |EAVerage|EPOSitive|MPOSitiv|RMS
%             fprintf(s.sa,'BAND 120kHz') %set RBW 120K
%             fprintf(s.sa,'BWID:VID:AUTO ON') % set VID auto.
%             fprintf(s.sa,':FREQ:STAR 0.1MHZ;\n') %just for now
%             fprintf(s.sa,':FREQ:STOP 2000MHZ;\n')  
%             fprintf(s.sa,'SWE:POIN 401\n')  
%             fprintf(s.sa,'POW:ATT 0 DB') % is undone by enabling preselector 8(
%             fprintf(s.sa,'INPut:PRESelector:PATH FILTer')
%             fprintf(s.sa,'POW:PRES:ATT 0 DB')
%             fprintf(s.sa,'POW:PRES:GAIN ON')
%             % %fprintf(sa,':MMEMory:LOAD:STAT 1,''C:\DTDFILES\STATE000.STA''; *WAI\n')
%             % fprintf(sa,':MMEMory:LOAD:CORR ANT, ''C:\DTDFILES\3142C3MH.ANT''; *WAI\n')
%             % fprintf(sa,':MMEMory:LOAD:CORR ANT, ''C:\DTDFILES\3142C3MV.ANT''; *WAI\n')
%             % fprintf(sa,':MMEMory:LOAD:CORR CABL, ''C:\DTDFILES\CHAMBER.CBL''; *WAI\n')
%             % fprintf(sa,':MMEMory:LOAD:LIM LLINE1, ''C:\data\''; *WAI\n')
%             
