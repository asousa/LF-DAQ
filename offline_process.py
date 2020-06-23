from cmd import Cmd
from sys import exit,argv,stdin
import sys, os
from getopt import getopt
import os.path
from optparse import OptionParser
from xml.dom.minidom import parse
import numpy as np

from utilities.DAQConfig import DAQConfig
# from Engine import Engine
from PostProcessorTree import PostProcessorTree
from PostProcessors.Narrowband import Narrowband
from PostProcessors.NarrowbandC import NarrowbandC

from datetime import datetime, timedelta

from utilities.DAQLogger import DAQLogger, DAQLogClient
from utilities.DAQConfig import DAQConfig
from utilities.loadMATdata import loadMATdata

import time
from __ver__ import __ver__
    
class Engine:

    def __init__(self, options, args, log_queue):

        self.stationSettings = options.settings.GetFirstSubTree("StationSettings")
        self.pptSettings = options.settings.GetFirstSubTree("PostProcessorTree")
       
        self.logger = DAQLogClient(log_queue, "ENG")
        self.PPT_logger= DAQLogClient(log_queue, "ENG.PPT")
        self.log_queue = log_queue
        self.logger.status('********** Initializing Engine **********')

        self.start()
        self.logger.status("Constructed PostProcessorTree.")


    def start(self):
        ''' Start the PostProcessorTree '''
        self.ppt = PostProcessorTree(self.pptSettings, self.stationSettings, self)
        # self.nbproc = Narrowband(options, self.PPT_logger)

    def stop(self):
        ''' Stop the PostProcessorTree '''
        self.logger.info("Shutting down PostProcessor")
        self.ppt.Stop()
        self.ppt = None
        self.logger.info('PPT process stopped')

    def get_queue_status(self):
        ''' Get the current length of the process queues '''

        qstats = []
        for sequence in self.ppt.ppt:
            for p in sequence:
                # try:
                if (len(p) > 0) and hasattr(p[0],'queue'):
                    qstats.append(p[0].queue.qsize())
                else:
                    logger.debug("Process has no queue to check")
                # except:
                #     logger.warning("issues with get_queue_status")
        return qstats


    def process_files(self, file_prefixes):
        '''
            Run prerecorded continuous files through the PostProcessorTree.
            Input: [file_prefixes]: A list of prefix files to load 
                   (e.g., strip the "_000.mat" and "_001.mat" suffixes off the end.)
                   We recreate the filenames in here based on the number of
                   channels in the PostProcessorTree.
        '''

        num_channels = self.ppt.GetNumSequences()

        for f_prefix in file_prefixes:

            # Start execution timer
            tstart = time.time()

            s_list = []
            for ind in xrange(num_channels):

                fname = f_prefix + "_%03d.mat"%ind
                logger.status("loading file %s"%fname)
                lm = loadMATdata(fname)
                s = lm.loadAllData()
                s_list.append(s)
        
                print('Data has dimensions', np.shape(s['data']))
                print('sample rate is ',s['Fs'])


            L = 1   #[seconds]
            fs = int(s['Fs'])
            logger.info("file length is %d seconds"%(len(s['data'])/s['Fs']))
            # t = np.arange(int(L*fs),dtype=np.float)/fs

            sec_in_file = len(s['data'])/s['Fs']
            n_strides = int(sec_in_file/L)

            stride_length = int(s['Fs']*L)

            logger.info("n_strides: %d"%n_strides)
            logger.info("stride_length: %d"%stride_length)

            lat = s['latitude']
            lon = s['longitude']
            alt = s['altitude']
            Q   = int(s['gps_quality'])
            timestamp = datetime(year=int(s['start_year']), month=int(s['start_month']), day=int(s['start_day']),
                hour=int(s['start_hour']), minute=int(s['start_minute']), second=int(s['start_second']))

            logger.info("time = %s"%timestamp)
            logger.info("data has shape %s"%np.shape(s['data']))

            for i in xrange(n_strides):
                        # for i in xrange(n_strides):

                # logger.info("stride %d"%i)
                left = i*stride_length
                right= (i+1)*stride_length
                # rawBuffer = s['data'][left:right].transpose().astype('float32')
                rawBuffer = [x['data'][left:right].transpose().astype('float32') for x in s_list]

        #        data = [[chNData], [dt,[lat,lon,alt],[quality]], sampleRate]
                data = [rawBuffer, [timestamp, [lat, lon, alt], [Q]], fs]

                # logger.info("Processing for timestamp %s." % data[1][0])
                self.ppt.Process(data)
                timestamp += timedelta(seconds=L)

            qstat = self.get_queue_status()
            logger.info("Queue status: %s"%(qstat))

            while any(np.array(qstat) > 0):
                qstat = self.get_queue_status()
                logger.timestamping("Waiting for PostProcessor Queues: depth =  %s          "%(qstat))
                time.sleep(1)

            # Lap
            tstop = time.time()
            run_ratio = sec_in_file/(tstop - tstart)
            logger.status('Finished file prefix %s (%0.1f faster than realtime)'%(f_prefix, run_ratio))


if __name__=="__main__":

    import multiprocessing
    multiprocessing.freeze_support()


    from optparse import OptionParser
    from xml.dom.minidom import parse
    
    if not os.path.isdir('log'):
        os.makedirs('log')

    # ========================
    # Setup options parser
    # ========================
    parser=OptionParser(version= __ver__,
                        description='VLF data acquisition program. ')
    parser.add_option("--settings",dest='settings_file', default='DaqSettings.xml', type="string",
                      help="Settings (.xml) filename.  Default: %default")
    parser.add_option('--autostart', dest="autostart", default=False, action="store_true",
                      help="Start daq acquisition without further input (default %default)")
    parser.add_option('--debug', dest="debug_on", default=False, action="store_true",
                      help="Start daq acquisition with debug mode turned on")
    parser.add_option('--in_dir', '--in','--input','--inp_dir', dest="inp_dir", default="Continuous",type="string",
                      help="The input directory of Continuous files to post-process. Default: %default")

    (options,args) = parser.parse_args()

    #print options
    
    # ========================
    # Config and Settings checking
    # ========================
    if not os.path.isfile(options.settings_file):
        #raise '%s does not exist' % options.settings #Error: exceptions must be old-style classes... JCC 04/02/12  
        sys.exit('Error: %s does not exist' % options.settings) # quit program if setting file doesn't exist
        
    print 'VLF DAQ Software v%s' % __ver__
    print 'Settings file: %s' % options.settings_file
    
    print "Parsing Settings file..."

    # Remove any Continuous MatFileWriter entries from the PostProcessorTree
    f = parse(options.settings_file)
    for el in f.getElementsByTagName("PostProcessor"):
        if el.getAttribute("module")  in "MatFileWriter":
            
            # Check if IsSynoptic == 0:
            if (el.getElementsByTagName('IsSynoptic')[0].childNodes[0].data) in "0":
                print('Removing Continuous MatFileWriter')
                parent = el.parentNode
                parent.removeChild(el)

    # Instantiate the configuration module
    options.settings = DAQConfig(f)

    # Start the loggers
    main_logger = DAQLogger(options)
    main_logger.start()    
    logger = DAQLogClient(main_logger.log_queue, "MAIN")
    logger.debug("Main logger ready.")


    inp_dir = options.inp_dir

    logger.status("Input directory:\t %s"%inp_dir)

    # ------- Engine version (using the whole PostProcessor tree) -----
    # # Start it!
    eng = Engine(options, args, main_logger.log_queue)


    if not os.path.isdir(inp_dir):
        logger.warning("Invalid input directory: %s"%inp_dir)
        engine.stop()
    else:

        # Find file prefixes to do:
        d = os.listdir(inp_dir)
        fnames = np.unique([x[:-8] for x in d if x.endswith('.mat')])
        prefixes = sorted([os.path.join(inp_dir, x) for x in fnames])
        
        # Run them!
        eng.process_files(prefixes)

        # Stop and start the engine, to finish the last .MAU files
        eng.stop()
        eng.start()
        eng.stop()

        logger.status("Finished processing files in %s"%inp_dir)