#!/bin/sh

#
# Get Raspi IP from internet if connected to the same network or from email if not
#
essid=`iwgetid -r`
if [ "$essid" = "Obi WAN" ]; then

    ip=`curl ifconfig.me 2> /dev/null`

else

    ip=`grep "Public IP" "$HOME/.thunderbird/e2l0inr0.default/Mail/Local Folders/Inbox.sbd/Raspi" | tail -n 1 | awk '{print $NF}'`

fi

#
# Replace old IP with new one
#
oldip=`grep -A 1 " pi" ~/.ssh/config | tail -1 | awk '{print $2}'`
line_no=`grep -nA 1 " pi" ~/.ssh/config | tail -1 | cut -d "-" -f1`

if [ "$ip" = "" ]; then
    exit
fi

if [ "$ip" != "$oldip" ]; then
  sed -i "${line_no}s/HostName.*$/HostName $ip/" ~/.ssh/config
fi
