# LF DAQ Distro -- CU LAIR Edition
#### Austin Sousa
#### 6.24.2018 - 6.5.2020
#### austin.sousa@colorado.edu


Hey! Here's a few notes about this version of the LF-DAQ distribution.

#### To run the DAQ:
- To run it from the bare python, use "run_VLF.bat" or "run_LF.bat"
- The .bat files automatically build the daqSettings xml file, using "Default_LF_settings.txt" or "Default_VLF_settings.txt"

#### To build a standalone DAQ:
- To build a standalone version, run "build_standalone.bat"
	- A new version of the distro will be packed up into the "Builds" directory.
	- You can copy the new build to a Windows machine without having to install Python or 
	  associated libraries.
	- You'll still need to install NIDAQMX - tested with version 17.1

#### Some new settings:
- copy_to_external:
	- if active=1, data will be stored in data_root on the internal drive
	- every 30 minutes, completed files will be moved to the external drive.
	- external_drive is a comma-separated list of external drives to use:
		external_drive=D:\,E:\,F:\
	- the DAQ software will cycle through this list, using the first drive until it's full,
		then moving on to the next, and so forth.

- copy_to_dropbox:
	- The intent of this module is to copy current spectrograms and the log file to dropbox,
		so it can be monitored remotely (without having to use the Stanford SSH server)
	- The transfer occurs every copy_period seconds
	- if tail_only=0, the entire log file will be copied to dropbox
	- if tail_only=1, only the last ~50 lines will be copied over, to save on bandwidth.
		(the log files can get huge)

- ErrorEmail:
	- Same as before -- every time the daq throws a major error, it'll send an email.
	- the email settings are brought out into the config file

- ErrorPost:
	- This *probably* works, but I haven't tested it since I can't get at the
	 old Stanford servers anymore.

#### To perform a calibration:
	- Disconnect the antenna loops from the preamp, so that we aren't recording any natural signals.
	- Adjust the configuration file to write 60-second files, and one 10-second spectrogram every minute.
	- Start the DAQ.
	- wait until you're taking data within a full-sized file (e.g., wait until the time is at :00:00)
	- Wait a few more seconds, then push the caltone button on the preamp box, and hold it for ~10 seconds.
	- quit the DAQ after the current file finishes. You should see the comb tone in the spectrograms.
	- run the calibration script: 
		- run "calibrate.bat", or the calibration .exe file in "software".
		- The script will ask you to open one of the data files (_000.mat, _001.mat).
		- The script will try to find the caltone in the data by looking for a large signal. You can also
		  manually tell the script where the caltone and a quiet region is, if the system is very noisy.
		- The script will ask you some details about the antenna. You can hit enter to use the defaults,
		  or enter new values for different antennas.
			- Square or triangle
			- wire AWG
			- number of turns
			- baseline, in centimeters.

			For RELAMPAGO, we're using:
				triangle
				16 AWG
				13 turns
				260 cm baseline

		- If it all went right, four plots will pop up, showing:
			- CalibrationResponse: 	The frequency response 
			- CalibrationNumber: 	The scaling factor (pT per amplitude division in sampling)
			- ResponseRatio: 		The matching ratio between the two channels
			- NoiseResponse: 		The frequency spectrum of the noise floor

		- The script will then ask you where to save the final files, and give you an option to
		  log any metadata (temporary calibration, new hardware, etc)
		- The computed calibration curves are stored in CalibrationVariables.mat
		- That's it! Close the plots, and copy CalibrationVariables.mat to the root directory.
		  The spectrogram module will automatically load it the next time the DAQ is run.

