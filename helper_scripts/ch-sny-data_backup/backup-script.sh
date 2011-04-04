#!/bin/bash

# This script does personal backups to a rsync backup server. You will end up
# with a 7 day rotating incremental backup. The incrementals will go
# into subdirectories named after the day of the week, and the current
# full backup goes into a directory called "current"
# Adaption of a script by tridge@linuxcare.com,
# made by Yun William Yu <yun.yu09@imperial.ac.uk>

# directories to backup; make sure there is NO trailing slash!
#BDIR="/home/$USER /media/STORE/lab"
BDIR="/home/$USER"

# excludes file - this contains a wildcard pattern per line of files to exclude
EXCLUDES=$HOME/cron/excludes

# the name of the backup machine
BSERVER=ch-sny-data.ch.ic.ac.uk

# Should be an absolute path to the backup directory on the server
BACKUPPATH=/home/$USER/backup

########################################################################

BACKUPDIR=`date +%A`
OPTS="--force --ignore-errors --delete-excluded --exclude-from=$EXCLUDES 
      --delete --backup --backup-dir=$BACKUPPATH/$BACKUPDIR -a "

export PATH=$PATH:/bin:/usr/bin:/usr/local/bin

# the following line clears the last weeks incremental directory
[ -d $HOME/cron/emptydir ] || mkdir $HOME/cron/emptydir
rsync --delete -a $HOME/cron/emptydir/ $BSERVER:$BACKUPPATH/$BACKUPDIR/ -e "ssh -i $HOME/.ssh/backup"
rmdir $HOME/cron/emptydir

# now the actual transfer
rsync $OPTS $BDIR $BSERVER:$BACKUPPATH/current -e "ssh -i $HOME/.ssh/backup" 

