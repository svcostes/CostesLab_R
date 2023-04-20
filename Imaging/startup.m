
fprintf(">>>LOADING startup.m for genuser\n")

fprintf(">>>LOADING DIPLIB CODEBASE\n")
addpath('/data/install_diplib/dip/common/dipimage') % version 2.9 import
fprintf(">>>INITIALIZING DIPLIB\n")
dip_initialise % also for version 2.9

image_codebase_dir = '/home/genuser/local_repos/Imaging';
fprintf(">>>LOADING IMAGING CODE AT %s \n", image_codebase_dir)
addpath(image_codebase_dir) % a set of scripts originally imported to this machine by Sylvain, used for image processing

addpath('/home/genuser/Documents/MATLAB/bfmatlab'); % path to have bio-format installed


% custom user settings
setenv("EDITOR", "sudo vim") % here to make the default editor when launching from CLI vim.  Without setting this, it tries to open the GUI editor which doesn't work if you are ssh'ed in.

fprintf(">>>DONE LOADING /home/genuser/startup.m\n")
fprintf(">>>There should be NO ERRORS ABOVE THIS LINE, OTHERWISE startup.m MAY HAVE ISSUES\n\n")

