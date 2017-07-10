#!/bin/sh

ip=`curl ifconfig.me 2> /dev/null`
oldip=`grep -A 1 " pi" ~/.ssh/config | tail -1 | awk '{print $2}'`

if [ "$ip" != "$oldip" ]; then
  sed -i "s/HostName.*$/HostName $ip/" ~/.ssh/config
fi
