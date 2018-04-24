"""
Copy files to a directory
"""
import shutil

__all__ = ['CopyFiles']

from Task import *

try:

    import wx
    from VLFPanel import VLFPanel


    class GUIPanel(VLFPanel):
        def __init__(self, parent):
            VLFPanel.__init__(self, parent, "CopyFiles")
            self.widgets = {
                'WriteDirectory':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'Interval':(wx.TextCtrl(self, wx.ID_ANY, size= self.FLOAT_SIZE),'[Seconds]'),
                'Delete':wx.Choice(self, wx.ID_ANY, size= self.DROP_BOOL_SIZE,choices = ['0','1']),
                'Filename':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'Filename1':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'Filename2':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'Filename3':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'Filename4':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'Swap':wx.TextCtrl(self, wx.ID_ANY, size= self.STRING_SIZE),
                'Directory':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'Directory1':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'Directory2':wx.TextCtrl(self, wx.ID_ANY, size= self.DIR_SIZE),
                'StartTime':(wx.TextCtrl(self, wx.ID_ANY, size= self.FLOAT_SIZE),'[HH:MM], System UTC'),
                'EndTime':(wx.TextCtrl(self, wx.ID_ANY, size= self.FLOAT_SIZE),'[HH:MM], System UTC'),
            }
            self.addWidgets()
except:
    print "WARNING: CAN'T IMPORT WX AND/OR VLFPanel"

class CopyFiles(Task):

###############################################################################

    def __init__(self, config, logger, tasknum):
    
        Task.__init__(self, config, logger, tasknum)

        #broadband, jpeg, and most narrowband files:
        self.datetime_re = re.compile('.*(\d{2})(\d{2})(\d{2})\d{6}[A-Za-z]{0,3}_\d{3}[A-Fa-f]{0,1}\..+')

        self.taskname = 'T%d.CopyFiles' % tasknum

        self.exclude_list = ['.mau']

        #parse the xml file
        self.writeDirectory = config.GetStrElemVal('WriteDirectory')
        if not os.path.isdir(self.writeDirectory):
            try:
                os.makedirs(self.writeDirectory)
            except:
                self.logger.error('Cannot create folder %s for copying--drive may be missing or full' % self.writeDirectory)

        self.interval = config.GetIntElemVal('Interval') #seconds
        self.delete_after_copy = config.GetIntElemVal("Delete",0)>0

        self.filenames = []
        no_file_default = 'X3@X-+#@!3blablabla'
        for ii in ['','0','1','2','3','4','5','6','7','8','9']:
            read_file = config.GetStrElemVal("Filename%s" % ii,no_file_default)
            if read_file != no_file_default:
                self.filenames.append(read_file)


        self.directories = []
        no_dir_default = 'X3@X-+#@!3blablabla'
        for ii in ['','0','1','2','3','4','5','6','7','8','9']:
            read_dir = config.GetStrElemVal("Directory%s" % ii,no_dir_default)
            if read_dir != no_dir_default:
                if os.path.isdir(read_dir):
                    self.directories.append(read_dir)
                else:
                    print 'Directory %s does not exist; skipping' % read_dir



###############################################################################

    def DoTask(self):

        #check write directory:
        if not os.path.isdir(self.writeDirectory):
            try:
                os.makedirs(self.writeDirectory)
            except:
                self.logger.error('Cannot create write directory %s' % self.writeDirectory)
                return

        self.CopyFiles(self.filenames)

        for directory in self.directories:
            if not self.running:
                break   #don't wait for everything to copy when software called to exit

            if not os.path.isdir(directory):
                continue

            #get a list of all files in the data directory
            file_list = os.listdir(directory)

            #prepend directory
            file_list = [os.path.join(directory,fname) for fname in file_list]

            #strip out all directories from the list
            file_list = filter(lambda x: not os.path.isdir(x), file_list)

            #strip out certain extensions
            for extension in self.exclude_list:
                file_list = filter(lambda x: not x.endswith(extension), file_list)

            try:
                self.CopyFiles(file_list)
                self.logger.info("Finished copying files.")
            except:
                # Would be nice to differentiate between different file errors (ie disc full, etc) and handle accordingly
                continue

    def CopyFiles(self,file_list):
        for filename in file_list:
            if not self.running:
                break   #don't wait for everything to copy when software called to exit
            try:
                if os.path.isfile(filename):
                    (head,tail) = os.path.split(filename)

                    #check for regular expression: date string in file
                    """
                    m = self.datetime_re.match(tail)
                    if m is not None:
                        new_dir = os.path.join(self.writeDirectory,"20%s" % m.group(1),"%s" % m.group(2),"%s" % m.group(3))
                        if not os.path.isdir(new_dir):
                            os.makedirs(new_dir)
                        new_filename = os.path.join(new_dir,tail)
                    else:
                        new_filename = os.path.join(self.writeDirectory,tail)
                    """
                    ### changed 4/21/13 by JCC to make copy destination directory flat
                    new_filename = os.path.join(self.writeDirectory,tail)

                    shutil.copyfile(filename,new_filename)
                    self.logger.debug("Copied %s to %s." % (filename,new_filename))
                    
                    if self.delete_after_copy:
                        os.remove(filename)
                        self.logger.debug("Removed %s." % filename) 
                        
            except Exception:
                self.logger.exception('Filename: %s' % filename)
                #break   #usually write drive dne; exit out of loop
                raise


###############################################################################

if __name__ == '__main__':
    from xml.dom.minidom import parseString
    import time as timemod
    import logging
    logger = logging.getLogger('')

    f = open('hi.txt','w')
    f.close()

    f = open('hi1.txt','w')
    f.write('phase 1')
    f.close()

    timemod.sleep(2)

    # Create the Motorola GPS clock object
    settings = """\
    <Task module="CopyFiles">
        <WriteDirectory>hi</WriteDirectory>
        <Interval>5</Interval>
        <Filename>hi.txt</Filename>
        <Filename1>hi1.txt</Filename1>
        <Filename2>hi2.txt</Filename2>
        <Directory>dir1</Directory>
        <Directory1>dir2</Directory1>


        <Delete>1</Delete>
    </Task>
    """
    copy = CopyFiles(parseString(settings),logger)
    copy.Start()
    timemod.sleep(8)
    f = open('hi1.txt','w')
    f.write('phase 2')
    f.close()
    timemod.sleep(6)
    copy.Stop()

