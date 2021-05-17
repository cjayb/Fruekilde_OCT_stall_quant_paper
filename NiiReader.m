classdef NiiReader < handle
    % NiiReader object for reading OCT files in NIFTI format
    %
    % Class instances wrap methos for opening and efficient reading of
    % large 3D OCT files saved in Nifti format. Basic usage:
    % >>> nii = NiiReader(fullfile(pathname, filename));
    % >>> frame = nii.readFrame(150);  % read 150th XZ-plane
    % >>> imagesc(frame)  % already converted to reflexivity
    properties
        fid  % file ID for open files
        dims  % dimensions of volume (x, y, z)
        info  % detailed information on volume contents
        zrange  % (optional) sub-range to read in the z-dimension
    end
    methods
        function obj = set.zrange(obj, value)
            obj.zrange = value;
        end
        function obj = NiiReader(in_file)
            % NiiReader constructor method
            %
            % Argument: input nii-file (string)
            % Returns: the reader object

            [fid, nk, nx, ny, info] = obj.open_nii(in_file);
            obj.fid = fid;
            obj.dims = [nk, nx, ny];
            obj.info = info;
            obj.zrange = false;
        end

        function obj = close(obj)
            % close the Nifti file after reading
            obj.fid = fclose(obj.fid);
            obj.dims = [];
            obj.info = [];
        end
        function ref_frame = readFrame(obj, yslice)
            % readFrame read a 2D slice from disk and convert spectrum to
            % reflectivity profile using NII2RR.m from Boas lab.
            %
            % Arguments: yslice specifies the index in the y-dim to read
            % Returns: the 2D slice read
            %
            % NB: if the zrange-property is set, only the limited z-range is
            % returned
            nk = obj.dims(1);
            nx = obj.dims(2);
            ny = obj.dims(3);
            if yslice > ny
                error('y-dimension is %d, you asked for slice %d', ny, yslice)
            end
            fseek(obj.fid,(yslice - 1) * nk * nx * 2 + 352, 'bof');
            frame_data = fread(obj.fid, nk * nx, 'int16');
            frame = reshape(frame_data, [nk nx]);
            
            ref_frame = NII2RR(single(frame));  % NB needed in path

            if obj.zrange
                ref_frame = ref_frame(obj.zrange(1):obj.zrange(2), :);
            end
        end
    end
   
    methods(Static)
        function [fid, nk, nx, ny, info] = open_nii(fname)
            
            [info] = NiiReader.nii_read_header(fname);
            % NB we only really need the dimensions, the rest of the header
            % information is not used. But reading the whole thing acts as a
            % form of unit test, so keeping it in!
            dim = info.Dimensions(1,1:3);
            nk = dim(1); nx = dim(2); ny = dim(3);
            
            % read data
            fid=fopen(fname,'rb');   
        end
        function [info] = nii_read_header(filename)
            % function for reading header of NifTi ( .nii ) volume file
            %
            % info = nii_read_header(filename);
            %
            % examples:
            % 1,  info=nii_read_header()
            % 2,  info=nii_read_header('volume.nii');

            if(exist('filename','var')==0)
                [filename, pathname] = uigetfile('*.nii', 'Read nii-file');
                filename = [pathname filename];
            end
            
            bswap=false;
            test=true;
            while(test)
                if(bswap)
                    fid=fopen(filename,'rb','b');
                else
                    fid=fopen(filename,'rb','l');
                end
                if(fid<0)
                    fprintf('could not open file %s\n',filename);
                    return
                end
                
                %get the file size
                fseek(fid,0,'eof');
                info.Filesize = ftell(fid);
                fseek(fid,0,'bof');
                info.Filename=filename;
                info.SizeofHdr=fread(fid,1,'int');
                info.DataType=fread(fid, 10, 'uint8=>char')';
                info.DbName=fread(fid, 18, 'uint8=>char')';
                info.Extents=fread(fid,1,'int');
                info.SessionError=fread(fid,1,'uint16');
                info.Regular=fread(fid, 1, 'uint8=>char')';
                info.DimInfo=fread(fid, 1, 'uint8=>char')';
                swaptemp=fread(fid, 1, 'uint16')';
                info.Dimensions=fread(fid,7,'uint16')'; % dim = [ number of dimensions x,y,z,t,c1,c2,c3];
                
                if(swaptemp(1)<1||swaptemp(1)>7), bswap=true; fclose(fid); else test=false; end
            end
            info.headerbswap=bswap;
            info.IntentP1=fread(fid,1,'float');
            info.IntentP2=fread(fid,1,'float');
            info.IntentP3=fread(fid,1,'float');
            info.IntentCode=fread(fid,1,'uint16');
            info.DataType=fread(fid,1,'uint16');
            datatypestr{1}={0,'UNKNOWN',  0}; % what it says, dude
            datatypestr{2}={1,'BINARY',   1}; % binary (1 bit/voxel)
            datatypestr{3}={2,'UINT8'  ,  8};% unsigned char (8 bits/voxel)
            datatypestr{4}={4,'INT16'   , 16}; % signed short (16 bits/voxel)
            datatypestr{5}={8,'INT32'  ,  32}; % signed int (32 bits/voxel)
            datatypestr{6}={16,'FLOAT' ,  32}; % float (32 bits/voxel)
            datatypestr{7}={32,'COMPLEX', 64}; % complex (64 bits/voxel)
            datatypestr{8}={64,'DOUBLE',  64}; % double (64 bits/voxel)
            datatypestr{9}={128,'RGB'  ,  24}; % RGB triple (24 bits/voxel)
            datatypestr{10}={255,'ALL'  ,  0}; % not very useful (?)
            datatypestr{11}={256,'INT8' ,  8}; % signed char (8 bits)
            datatypestr{12}={512,'UINT16', 16}; % unsigned short (16 bits)
            datatypestr{13}={768,'UINT32', 32}; % unsigned int (32 bits)
            datatypestr{14}={1024,'INT64', 64}; % long long (64 bits)
            datatypestr{15}={1280,'UINT64',     64}; % unsigned long long (64 bits)
            datatypestr{16}={1536,'FLOAT128',   128}; % long double (128 bits)
            datatypestr{17}={1792,'COMPLEX128', 128}; % double pair (128 bits)
            datatypestr{18}={2048,'COMPLEX256', 256}; % long double pair (256 bits)
            datatypestr{19}={2304,'RGBA32', 32}; % 4 byte RGBA (32 bits/voxel)
            info.datatypestr='UNKNOWN';
            info.bitvoxel=0;
            for i=1:19
                if(datatypestr{i}{1}==info.DataType)
                    info.DataTypeStr=datatypestr{i}{2};
                    info.BitVoxel=datatypestr{i}{3};
                end
            end
            
            info.Bitpix=fread(fid,1,'uint16');
            info.SliceStart=fread(fid,1,'uint16');
            temp=fread(fid,1,'float');
            info.PixelDimensions=fread(fid,7,'float');
            
            info.VoxOffset=fread(fid,1,'float');
            info.RescaleSlope=fread(fid,1,'float');
            info.RescaleIntercept=fread(fid,1,'float');
            info.SliceEnd=fread(fid,1,'uint16');
            info.SliceCode=fread(fid, 1, 'uint8=>char')';
            info.XyztUnits=fread(fid, 1, 'uint8')';
            dataunitsstr{1}={'UNKNOWN', 0}; %! NIFTI code for unspecified units.
            dataunitsstr{2}={'METER',   1};  %! NIFTI code for meters.
            dataunitsstr{3}={'MM',    2};  %! NIFTI code for millimeters.
            dataunitsstr{4}={'MICRON ', 3};  %! NIFTI code for micrometers.
            dataunitsstr{5}={'SEC',    8};  %! NIFTI code for seconds.
            dataunitsstr{6}={'MSEC',   16};  %! NIFTI code for milliseconds.
            dataunitsstr{7}={'USEC',  24};  %! NIFTI code for microseconds.
            dataunitsstr{8}={'HZ',  32};  %! NIFTI code for Hertz.
            dataunitsstr{9}={'PPM',  40};  %! NIFTI code for ppm.
            dataunitsstr{10}={'RADS',  48};  %! NIFTI code for radians per second.
            info.xyzt_unitsstr='UNKNOWN';
            for i=1:10,
                if(dataunitsstr{i}{2}==info.XyztUnits)
                    info.XyztUnitsStr=dataunitsstr{i}{1};
                end
            end
            
            info.CalMax=fread(fid,1,'float');
            info.CalMin=fread(fid,1,'float');
            info.Slice_duration=fread(fid,1,'float');
            info.Toffset=fread(fid,1,'float');
            info.Glmax=fread(fid,1,'int');
            info.Glmin=fread(fid,1,'int');
            info.Descrip=fread(fid, 80, 'uint8=>char')';
            info.AuxFile=fread(fid, 24, 'uint8=>char')';
            info.QformCode=fread(fid,1,'uint16');
            info.SformCode=fread(fid,1,'uint16');
            info.QuaternB=fread(fid,1,'float');
            info.QuaternC=fread(fid,1,'float');
            info.QuaternD=fread(fid,1,'float');
            info.QoffsetX=fread(fid,1,'float');
            info.QoffsetY=fread(fid,1,'float');
            info.QoffsetZ=fread(fid,1,'float');
            info.SrowX=fread(fid,4,'float');
            info.SrowY=fread(fid,4,'float');
            info.SrowZ=fread(fid,4,'float');
            info.IntentName=fread(fid, 16, 'uint8=>char')';
            info.Magic=fread(fid, 4, 'uint8=>char')';
            
            fclose(fid);
        end
    end
end