#!/bin/sh

ip=`curl ifconfig.me 2> /dev/null`
oldip=`grep -A 1 " pi" ~/.ssh/config | tail -1 | awk '{print $2}'`
line_no=`grep -nA 1 " pi" ~/.ssh/config | tail -1 | cut -d "-" -f1`

if [ "$ip" != "$oldip" ]; then
  sed -i "${line_no}s/HostName.*$/HostName $ip/" ~/.ssh/config
fi
