classdef edfinfo < matlab.mixin.CustomDisplay & ...
        matlab.mixin.SetGet 
%edfinfo Get information from the EDF or EDF+ file header
%   INFO = EDFINFO(FILENAME) returns header information of the
%   corresponding European Data Format (EDF or EDF+) file with location
%   specified in string FILENAME.
%   INFO is an edfinfo object with these read-only properties:
%   'Filename'           - Name of the file.
%   'FileModDate'        - A string that specifies the modification date
%                          of the file.
%   'FileSize'           - An integer indicating the size of the file in
%                          bytes.
%   'Version'            - A string that specifies the data format version
%                          and whose value is always "0".
%   'Patient'            - A string that contains patient identification
%                          details such as:
%                               - Patient ID
%                               - Gender
%                               - Birth Date
%                               - Patient Name
%   'Recording'          -  A string that contains recording identification
%                           details such as:
%                               - Start date/time
%                               - Technician ID
%                               - Equipment ID
%   'StartDate'          - A string that specifies the start date of the
%                          recording in 'dd.MM.yy' format.
%   'StartTime'          - A string that specifies the start time of the
%                          recording in 'HH.mm.ss' format.
%   'HeaderBytes'        - An integer indicating the size of the header in
%                          bytes.
%   'Reserved'           - A string that indicates if the file is
%                          continuous ("EDF+C") or discontinuous ("EDF+D"),
%                          when the file is EDF+ compliant. If the file
%                          type is not EDF+, then 'Reserved' is empty.
%   'NumDataRecords'     - An integer indicating the number of data records
%                          present in the specified file.
%   'DataRecordDuration' - A duration scalar indicating the duration of
%                          each data record in seconds.
%   'NumSignals'         - An integer indicating the number of signals
%                          present in the specified file.
%   'SignalLabels'       - A 'NumSignals'-by-1 string vector of signal
%                          names.
%   'TransducerTypes'    - A string vector containing details about the
%                          transducers used to record the data. Each
%                          element of 'TransducerTypes' corresponds to a
%                          signal in 'SignalLabels'.
%   'PhysicalDimensions' - A string vector containing the physical
%                          dimensions or units of the recorded data. Each
%                          element of 'PhysicalDimensions' corresponds to a
%                          signal in 'SignalLabels'.
%   'PhysicalMin'        - A 'NumSignals'-by-1 array indicating the minimum
%                          physical value of each recorded signal.
%   'PhysicalMax'        - A 'NumSignals'-by-1 array indicating the maximum
%                          physical value of each recorded signal.
%   'DigitalMin'         - A 'NumSignals'-by-1 array indicating the minimum
%                          digital value of each recorded signal.
%   'DigitalMax'         - A 'NumSignals'-by-1 array indicating the maximum
%                          digital value of each recorded signal.
%   'Prefilter'          - A string vector containing details of the type
%                          of filter applied when recording the signals.
%                          Each element of 'Prefilter' corresponds to a
%                          signal in 'SignalLabels'.
%   'NumSamples'         - A 'NumSignals'-by-1 array indicating the number
%                          of samples of the signals present in each
%                          data record.
%   'SignalReserved'     - A 'NumSignals'-by-1 string vector containing
%                          any additional information about the signals.
%   'Annotations'        - A timetable containing all annotations present
%                          in the data records. Each row represents an
%                          annotation in the file. The timetable contains
%                          these variables:
%
%                        Onset      - Time at which the annotation
%                                     started, specified as a duration
%                                     indicating the number of seconds
%                                     relative to the start time of the
%                                     file.
%                        Annotation - A string containing the annotation
%                                     text.
%                        Duration   - A duration scalar indicating the
%                                     duration of the event described
%                                     by the annotation. If the file does
%                                     not specify an annotation duration,
%                                     this variable is returned as
%                                     NaN.
%
%                        "Annotations" is an empty timetable when there
%                         are no annotations in the EDF or EDF+ file.
%
%   % EXAMPLE:
%      % Read header of example.edf
%      info = edfinfo('example.edf')
%
%   See also EDFREAD, EDFHEADER, EDFWRITE

%   Copyright 2020 The MathWorks, Inc.

%   References:
% 	  [1] Bob Kemp, Alpo Värri, Agostinho C. Rosa, Kim D. Nielsen, and
%         John Gade. "A simple format for exchange of digitized polygraphic
%         recordings." Electroencephalography and Clinical
%         Neurophysiology 82 (1992): 391-393.
% 	  [2] Bob Kemp and Jesus Olivan. "European data format 'plus' (EDF+),
%         an EDF alike standard format for the exchange of physiological
%         data." Clinical Neurophysiology 114 (2003): 1755-1761.
    
properties (SetAccess = private)
    Filename
    FileModDate
    FileSize
    Version
    Patient
    Recording
    StartDate
    StartTime
    HeaderBytes
    Reserved
    NumDataRecords
    DataRecordDuration
    NumSignals
    SignalLabels
    TransducerTypes
    PhysicalDimensions
    PhysicalMin
    PhysicalMax
    DigitalMin
    DigitalMax
    Prefilter
    NumSamples
    SignalReserved
    Annotations
end

properties (Access = private)
    % Property list used to display the object
    pPropertyList = [...
        "Filename"...
        "FileModDate"...
        "FileSize"...
        "Version"...
        "Patient"...
        "Recording"...
        "StartDate"...
        "StartTime"...
        "HeaderBytes"...
        "Reserved"...
        "NumDataRecords"...
        "DataRecordDuration"...
        "NumSignals"...
        "SignalLabels"...
        "TransducerTypes"...
        "PhysicalDimensions"...
        "PhysicalMin"...
        "PhysicalMax"...
        "DigitalMin"...
        "DigitalMax"...
        "Prefilter"...
        "NumSamples"...
        "SignalReserved"...
        "Annotations"];
end

methods (Hidden)
    function obj = edfinfo(filename)
        % Open the EDF/EDF+ file
        [obj, filename, fid, fileInfo] = openFile(obj, filename);
        
        % Close the opened file using onCleanup
        cleanup = onCleanup(@() fclose(fid));
        
        % Read header of EDF/EDF+ file
        obj = readHeader(obj, filename, fid, fileInfo);
    end
end

%----------------------------------------------------------------------
% readHeader methods
%----------------------------------------------------------------------
methods (Access = private)

    function obj = assignValues(obj,props,values)
        for ii = 1:numel(props)
            obj.(props{ii}) = values{ii};
        end
    end

    function [obj,filename,fid,fileInfo] = openFile(obj,filename)
        % Convert string arrays to character arrays
        filename = convertStringsToChars(filename);
        
        validateattributes(filename,{'char','string'},...
            {'scalartext','nonempty'},mfilename);
        
        % Error out when the file extension is not .edf/.EDF
        [~, ~, ext] = fileparts(filename);
        if ~strcmpi(ext,'.edf')
            error(message('signal:edf:InvalidFileExt'));
        end
        
        % Open the EDF file
        [fid, fileInfo] = signal.internal.edf.openFile(filename,'r');
    end
    
    function obj = readHeader(obj,filename,fid,fileInfo)
        
       [version, patient, recording, startDate, startTime, headerBytes,...
        reserve, numDataRecords, dataRecordDuration, numSignals,...
        sigLabels, transducerType, physicalDimension, physicalMinimum,...
        physicalMaximum, digitalMinimum, digitalMaximum, prefilter,...
        numSamples, sigReserve, ~, annotations,...
        ~, ~] = signal.internal.edf.edfinfo('edfinfo', fid, ...
        filename, fileInfo, false);
        
        % Project matching properties from internal info onto this obj
        props = properties('edfinfo');
        
        values = {string(fileInfo.name), string(datestr(fileInfo.datenum)),... 
            fileInfo.bytes, string(version), string(patient),...
            string(recording), string(startDate), string(startTime),...
            headerBytes, reserve, numDataRecords, ...
            seconds(dataRecordDuration), numSignals, sigLabels,...
            transducerType, physicalDimension, physicalMinimum,...
            physicalMaximum, digitalMinimum, digitalMaximum, prefilter,...
            numSamples, sigReserve, annotations};
        
        obj = assignValues(obj, props, values);
    end
end

%----------------------------------------------------------------------
% Display control, Copy
%----------------------------------------------------------------------
methods (Access = protected)
    function propgrp = getPropertyGroups(obj)
        %getPropertyGroups Group properties in order for object display
        propList = obj.pPropertyList;
        propgrp = matlab.mixin.util.PropertyGroup(propList);
    end
end % end methods

end
