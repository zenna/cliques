#!/bin/bash
# install.sh
# Script to set up nightly backup of home directory

server=ch-sny-data.ch.ic.ac.uk
ssh-keygen -t dsa -f $HOME/.ssh/backup -N ''

scp $HOME/.ssh/backup.pub $USER@$server:/home/$USER/
ssh $USER@$server "if [ ! -d backup ]; then mkdir backup; fi ; if [ ! -d .ssh ]; then mkdir .ssh ; chmod 700 .ssh ; fi ; mv backup.pub .ssh/ ; cd .ssh/ ; if [ ! -f authorized_keys ] ; then touch authorized_keys ] ; chmod 600 authorized_keys ; fi ; cat backup.pub >> authorized_keys ; rm backup.pub;"
# At this point, you'll be prompted to login with your current password (twice)

if [ ! -d $HOME/cron ];
then
	mkdir $HOME/cron/
fi

# Copy backup script
cp backup-script.sh $HOME/cron/
chmod 755 $HOME/cron/backup-script.sh

# Create rsync excludes file with network directory "homedir"
echo $HOME/homedir/ >> $HOME/cron/excludes

echo You need to add the following line to the crontab, which you will edit in the next dialogue:
echo "0 3 * * * $HOME/cron/backup-script.sh >/dev/null 2>&1"
read -p "Press enter to edit the crontab..."

crontab -e
# You need to add the following line (without the comment #, of course)
# 0 3 * * * /home/<username>/cron/backup-script.sh >/dev/null 2>&1
# This should change, of course, if your home directory is located somewhere
# different.

echo Depending on your system, you may need to logout/restart in order for ssh-agent to properly cache the new keyfile.
