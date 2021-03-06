"""
MatFileWriter creates MATLAB mat files and appends data to them.
"""

__all__ = ['MatFileWriter']
from struct import pack
from sys import byteorder
from datetime import datetime, timedelta
import os, os.path
from scipy.io.numpyio import fwrite, fread
import subprocess

# Constants
types = {"float64": [0,'=d'], "float32": [1,'=f'], "int32": [2,'=i'], "int16": [3,'=h'], "uint16": [4,'=H'], "uint8": [5,'=B']}


try:
	import wx
	from VLFPanel import VLFPanel

	from main import __ver__

	class GUIPanel(VLFPanel):
	    def __init__(self, parent):
	        VLFPanel.__init__(self, parent, "MatfileWriter")
	        self.widgets = {
	            'adc_channel_number':(wx.TextCtrl(self, wx.ID_ANY, size= self.FLOAT_SIZE),'0-999'),
	            'DirectoryRoot':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
	            'DirectoryRoot1':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
	            'DirectoryRoot2':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
	            'DirectoryRoot3':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
	            'Duration':(wx.TextCtrl(self, wx.ID_ANY, size= self.FLOAT_SIZE),'[Seconds]'),
	##            'IsSynoptic':(wx.TextCtrl(self, wx.ID_ANY, size= self.BOOL_SIZE),'Default: 0'),
	            'IsSynoptic':(wx.Choice(self, wx.ID_ANY, size= self.DROP_BOOL_SIZE,choices = ['0','1']),''),
	            'SynopticPeriod':(wx.TextCtrl(self, wx.ID_ANY, size= self.FLOAT_SIZE),'[Seconds]'),
	            'cal_factor':wx.TextCtrl(self, wx.ID_ANY, size= self.FLOAT_SIZE),
	            'system_call':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE)
	        }
	        self.addWidgets()
except:
	__ver__ = 'TEST'

class dummyLogger():
    def __init__(self,name='MFW'):
    	self.name = name
    def debug(self,a=''):
        print '%s: DEBUG: %s' % (self.name,a)
    def info(self,a=''):
        print '%s: INFO: %s' % (self.name,a)
    def warning(self,a=''):
        print '%s: WARNING: %s' % (self.name,a)
    def error(self,a=''):
        print '%s: ERROR: %s' % (self.name,a)
    def critical(self,a=''):
        print '%s: CRITICAL: %s' % (self.name,a)

class nullLogger():
    def __init__(self):
    	pass
    def debug(self,a=''):
        pass
    def info(self,a=''):
        pass
    def warning(self,a=''):
        pass
    def error(self,a=''):
        pass
    def critical(self,a=''):
        pass


class MatFileWriter:
	"""
	The MatFileWriter class is a Post Processor class that writes data to file
	in the format of a Level-4 MAT File.
	"""

	# ========================
	# Constructors/Destructors
	# ========================
	def __init__(self, config,logger=None,narrowband=None):
		"""
		Constructs the MatFileWriter object according to the supplied XML
		description.
		"""
		# Initialize data members
		self.filename = None
		self.file_renamed = True
		self.fileStart = None
		self.datarowsptr = None
		self.datarows = None
		self.MOPToffset = 0
		self.startTime = None
		if byteorder == 'big':
			self.MOPToffset = 1000
		if logger is None:
			self.logger = nullLogger()
		else:
			self.logger = logger
		self.writing_file = 0
		self.narrowband = None
		self.forceStartNewFile = None   #new



		#check write directories:
		self.directoryRoots = []
		self.overwrite_names = []
		no_dir_default = 'X3@X-+#@!3blablabla'
		for ii in ['','0','1','2','3','4','5','6','7','8','9']:
			write_dir = self.GetStrElemVal(config, "DirectoryRoot%s" % ii,no_dir_default)
			if write_dir != no_dir_default:
				(adir,overwrite) = self.parseDirectory(write_dir)
				if not os.path.isdir(adir):
					os.makedirs(adir)
				else:
					self.logger.info('Checking %s for .mau files ...' % adir)
					file_list = os.listdir(adir)
					file_list = filter(lambda x: x.endswith('.mau'), file_list)
					for afile in file_list:
						newname = os.path.join(adir,afile.replace('.mau','.mat'))
						if os.path.isfile(newname): os.remove(newname)
						os.rename(os.path.join(adir,afile),newname)
					if len(file_list)>0:
						self.logger.info('Found and renamed %d .mau files' % len(file_list))


				self.directoryRoots.append(adir)
				self.overwrite_names.append(overwrite)


		self.is_synoptic = self.GetIntElemVal(config, "IsSynoptic",0)
		self.synoptic_duration = None
		self.file_duration = 60

		self.synoptic_duration = self.GetIntElemVal(config, "Duration")
		if(self.is_synoptic):
			self.file_duration = self.GetIntElemVal(config, "SynopticPeriod")
		else:
			#during continuous recording, the synoptic period and recording duration
			#are the same
			self.file_duration = self.synoptic_duration


		# Read channel header information
		self.station_id             = self.GetStrElemVal(config, "station_id")
		self.adc_channel_number 	= self.GetIntElemVal(config, "adc_channel_number")
		self.adc_sn 				= self.GetStrElemVal(config, "adc_sn",'.')
		self.adc_type 				= self.GetStrElemVal(config, "adc_type",'.')
		self.antenna_bearings 		= self.GetStrElemVal(config, "antenna_bearings",'.')
		self.antenna_description	= self.GetStrElemVal(config, "antenna_description",'.')
		self.cal_factor 			= self.GetDblElemVal(config, "cal_factor",1.0)
		self.computer_sn 			= self.GetStrElemVal(config, "computer_sn",'.')
		self.gps_sn 				= self.GetStrElemVal(config, "gps_sn",'.')
		self.hardware_description 	= self.GetStrElemVal(config, "hardware_description",'.')
		self.is_broadband 			= self.GetIntElemVal(config, "is_broadband",1)
		self.station_description 	= self.GetStrElemVal(config, "station_description",'.')
		self.station_name 			= self.GetStrElemVal(config, "station_name",'X')
##		self.Version 				= self.GetStrElemVal(config, "VERSION",'.')
		self.Version = __ver__

		self.sys_call               = self.GetStrElemVal(config, "system_call",'.')


		if narrowband is not None:
			self.narrowband = 1
			self.is_broadband = 0
			self.is_amp = narrowband[0] #self.GetIntElemVal(narrowband, "is_amp")
			self.is_msk = narrowband[1] #self.GetIntElemVal(narrowband, "is_msk")
			self.T = narrowband[2]  #self.GetStrElemVal(narrowband, "T", 'T')
			self.filter_taps = narrowband[3]
			self.Fc = narrowband[4] #new
			self.call_sign = narrowband[5]




    # =======
    # Methods
    # =======
	def Process(self, data,dtype='int16'):
		"""
		Writes the received data to file as a Level-4 MAT File.
		"""

		if self.TimeToStartNewFile(data[1][0]) or self.forceStartNewFile:
			if self.forceStartNewFile:
				self.forceStartNewFile = False  #RKS
			self.rename_mau()   #change file extension of last file
			self.filenames = []
			for adir in self.directoryRoots:
				if self.narrowband is not None:
					self.filenames.append("%s/%s%s%s_%03d%s.mau" % \
						(adir, self.station_id, data[1][0].strftime("%y%m%d%H%M%S"),self.call_sign, self.adc_channel_number, self.T))
				else:
					self.filenames.append("%s/%s%s_%03d.mau" % (adir, self.station_id, data[1][0].strftime("%y%m%d%H%M%S"),self.adc_channel_number))
			self.file_renamed = False
			self.fileStart = data[1][0]
			self.WriteHeader(data)
			self.WriteDataHeader(dtype)

		#NOTE: On first acquisition period, synoptic recording may extend beyond normal cutoff time
		if(data[1][0]-self.fileStart < timedelta(seconds=self.synoptic_duration)):
			self.AppendData(data[0])
		else:   #rename synotpic
			self.rename_mau()
		return [self.filenames[0], data[1][0], self.file_duration, data[2]]


	def StartNewFile(self):
		self.forceStartNewFile = True


	# ==============
	# Helper Methods
	# ==============

	def WriteHeader(self, data):
		self.WriteVector('start_year', data[1][0].year, 'float64')
		self.WriteVector('start_month', data[1][0].month, 'float64')
		self.WriteVector('start_day', data[1][0].day, 'float64')
		self.WriteVector('start_hour', data[1][0].hour, 'float64')
		self.WriteVector('start_minute', data[1][0].minute, 'float64')
		self.WriteVector('start_second', data[1][0].second, 'float64')
		self.WriteVector('latitude', data[1][1][0], 'float64')
		self.WriteVector('longitude', data[1][1][1], 'float64')
		self.WriteVector('altitude', data[1][1][2], 'float64')
		self.WriteVector('Fs', data[2], 'float64')
		self.WriteVector('gps_quality',self.decode_gps_quality(data[1][2][0]),'uint8')
		self.WriteVector('adc_channel_number', self.adc_channel_number, 'float64')
		self.WriteVector('adc_sn', self.adc_sn, 'uint8')
		self.WriteVector('adc_type', self.adc_type, 'uint8')
		self.WriteVector('antenna_bearings', self.antenna_bearings, 'uint8')
		self.WriteVector('antenna_description', self.antenna_description, 'uint8')
		self.WriteVector('cal_factor', self.cal_factor, 'float64')
		self.WriteVector('computer_sn', self.computer_sn, 'uint8')
		self.WriteVector('gps_sn', self.gps_sn, 'uint8')
		self.WriteVector('hardware_description', self.hardware_description, 'uint8')
		self.WriteVector('is_broadband', self.is_broadband, 'float64')
		self.WriteVector('station_description', self.station_description, 'uint8')
		self.WriteVector('station_name', self.station_name, 'uint8')
		self.WriteVector('VERSION', self.Version, 'uint8')
		if self.narrowband is not None:
			self.WriteVector('is_amp', self.is_amp, 'float64')
			self.WriteVector('is_msk', self.is_msk, 'float64')
			self.WriteVector('Fc', self.Fc, 'float64')
			self.WriteVector('call_sign', self.call_sign, 'uint8')
			self.WriteVector('filter_taps', self.filter_taps, 'float64')

	def decode_gps_quality(self,quality):
		#The TrueTimeClock returns 50 or 51 depending on locked or unlocked
		# because the TT gps takes a long time to lock.  Decode here:
		if quality < 50:    #motorolla
			return '%d' % quality
		if quality == 50:   #truetime
			return 'UNLOCKED'
		if quality == 51:   #truetime
			return 'LOCKED'
		raise 'Unrecognized quality factor'

	def WriteVector(self, name, data, dtype, text=False):
		#If data is a string, need to call ord on each element:
		if dtype is 'uint8':
			data_convert = lambda x: ord(x)
		else:   # do nothing
			data_convert = lambda x: x

		for filename in self.filenames:
			f = open(filename, 'a+b')
			MOPT = self.MOPToffset + types[dtype][0]*10
			if text: MOPT += 1  								#?????????
			f.write(pack('=l', MOPT))							# MOPT
			try:	f.write(pack('=l', len(data)))				# Number of rows
			except:	f.write(pack('=l', 1))
			f.write(pack('=l', 1))								# Number of columns
			f.write(pack('=l', 0))								# Imaginary, should be 0
			f.write(pack('=l', len(name)+1))						# Length of name+1
			f.write(pack("=%is" % len(name), name) + chr(0))		# Entry name ending with a \0 character
			try:    lenData = len(data)
			except: lenData = 1
			if lenData > 1:
				for i in range(lenData):
					f.write(pack(types[dtype][1], data_convert(data[i])))
			else:
				f.write(pack(types[dtype][1], data_convert(data)))


			f.close()

	def WriteDataHeader(self,dtype='int16'):
		for filename in self.filenames:
			f = open(filename, 'a+b')
			name = "data"										# Name of the matfile entry
			MOPT = self.MOPToffset + types[dtype][0]*10
			f.write(pack('=l',MOPT))								# MOPT
			self.datarowsptr = f.tell()							# Save a pointer to the number of rows
			self.datarows = 0                                   # Save the number of rows in the data
			f.write(pack('=l', 0))								# Number of rows
			f.write(pack('=l', 1))								# Number of columns
			f.write(pack('=l', 0))								# Imaginary, should be 0
			f.write(pack('=l', len("data")+1))					# Length of name+1
			f.write(pack("=%is" % len("data"), "data") + chr(0))	# Entry name ending with a \0 character
			f.close()

	def AppendData(self, data):
		self.datarows += len(data)
		for filename in self.filenames:
			f = open(filename, 'r+b')
			f.seek(self.datarowsptr, 0)
			f.write(pack('=l', self.datarows))
			f.seek(0, 2)    #move pointer to end of file
			fwrite(f, data.size, data) #write a second's worth of data to disk
			f.close()

	def TimeToStartNewFile(self, timestamp):
		if self.fileStart is None:
			return True
		beginningOfDay = datetime(year=timestamp.year, month=timestamp.month, day=timestamp.day)
		if ((timestamp-beginningOfDay).seconds % self.file_duration) == 0:
			return True
		if (timestamp - self.fileStart) >= timedelta(seconds=self.file_duration):
			return True
		return False

	def GetIntElemVal(self, config, tag,default = None):
		if default is not None:
			if len(config.getElementsByTagName(tag))==0:
				return default
		child = config.getElementsByTagName(tag)[0].firstChild
		if child is None:
			return 0
		else:
##			print 'converting an int'
			return int(child.data)

	def GetStrElemVal(self, config, tag,default=None):
		if default is not None:
			if len(config.getElementsByTagName(tag))==0:
				return default
		child = config.getElementsByTagName(tag)[0].firstChild
		if child is None:
			return '.'
		else:
##			print 'converting an str'
			return str(child.data)

	def GetDblElemVal(self, config, tag,default=None):
		if default is not None:
			if len(config.getElementsByTagName(tag))==0:
				return default
		child = config.getElementsByTagName(tag)[0].firstChild
		if child is None:
			return 0.0
		else:
##			print 'converting an float'
			return float(child.data)

	def parseDirectory(self,dir_string):
		"""
		For each directory, there is an option to include a filename:
		<DirectoryRootX>SomeFolder,overwriteFilename</DirectoryRootX>
		If given, then the file name of the file written to this folder is overwritten
		with overwriteFilename.  If the usual filename (with the time stamp) is desired, then
		the format is:
		<DirectoryRootX>SomeFolder</DirectoryRootX>
		where X is either an empty string or a numeric (0--9)

		This method returns (directory, overwrite_filename)
		where overwrite_filename is None if no overwrite name specified
		"""
		index_comma = dir_string.rfind(',')
		if index_comma==-1: #just a directory
			return(dir_string,None)
		else:
			return (dir_string[:index_comma],dir_string[index_comma+1:])

	def rename_mau(self):
##		print self.overwrite_names

		if self.file_renamed is False:

			for i,afile in enumerate(self.filenames):
				if self.overwrite_names[i] is not None:
					(head,tail) = os.path.split(afile)
					filename = os.path.join(head,self.overwrite_names[i])
				else: #name to .mat
					filename = list(afile)
					filename[-1] = 't'
					filename = "".join(filename)
				if os.path.isfile(filename):
					os.remove(filename)
##				print afile, filename
				os.rename(afile,filename)
				self.filenames[i] = filename    #replace filename

			#call system call.  Pass in first filename
			self.sys_call_after_rename(self.filenames[0])

			self.file_renamed = True

	def sys_call_after_rename(self,filename):
		"""
		If pass in a string in sys_call, then this is called
		Place an @ to insert the filename (using the first DirectoryRoot)
		"""
		if self.sys_call == '.':
			return

		filename_index =  self.sys_call.rfind('@')
		if filename_index == -1:
			sys_call = self.sys_call
		else:
			sys_call = list(self.sys_call)
			sys_call.pop(filename_index)
			sys_call.insert(filename_index,filename)
			sys_call = "".join(sys_call)


##		print 'passing to the os: %s' % sys_call
		run = subprocess.Popen(sys_call,stdout=subprocess.PIPE)



# =========
# Unit Test
# =========
if __name__ == "__main__":
	from xml.dom.minidom import parseString
	import numpy
	from datetime import datetime, timedelta
	import time as timemod

	# Create the Motorola GPS clock object
	settings = """\
	<PostProcessor module="MatFileWriter">
        <station_id>Test</station_id>
		<DirectoryRoot>.</DirectoryRoot>
		<adc_channel_number>0</adc_channel_number>
		<system_call>Spectrogram\\CreateSpecgramJPG.exe @ .\\latestNS.mat 1 100</system_call>
		<Duration>10</Duration>
	</PostProcessor>
	"""
	settings = """\
	<PostProcessor module="MatFileWriter">
        <station_id>Test</station_id>
		<DirectoryRoot>.</DirectoryRoot>
		<adc_channel_number>0</adc_channel_number>
		<Duration>10</Duration>
	</PostProcessor>
	"""

	mfw = MatFileWriter(parseString(settings))

	sampleRate = 100000.0
	time = datetime(year=2008, month=6, day=28, hour=11, minute=59, second=50)
	for i in range(25):
##		rawBuffer = numpy.arange(sampleRate, dtype=numpy.int16) + numpy.random.random(sampleRate)
		rawBuffer = 100*numpy.sin(2.0*numpy.pi*numpy.arange(sampleRate)/sampleRate*10000.0)
		rawBuffer = rawBuffer.astype(numpy.int16)
		data = [rawBuffer, [time, [0, 0, 0], [10]],sampleRate]
		print "Processing for timestamp %s." % data[1][0]
		mfw.Process(data)
		time += timedelta(seconds=1)
##		timemod.sleep(1)
