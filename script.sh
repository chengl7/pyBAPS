#!/bin/bash

domain=`hostname --domain`
server="$2.${domain}"

N=$1
echo " "
#echo "server: ${server}"

#echo "input parameters: ${@:3}"

host=$(hostname --long)
echo "in test ${host}"

if [ "${host}" == "${server}" ]; then
    echo "server ${host}"
    python server.py $N ${server} "${@:3}"
else
    sleep 5
    echo "client ${host}"
    python worker.py $N ${server}
fi
echo " "
