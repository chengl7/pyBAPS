#!/bin/bash

domain=`hostname --domain`
server="$1.${domain}"

echo "server: ${server}"

host=$(hostname)
echo "in test $(hostname)"

if [ "$(hostname)" == "$1" ]; then
    echo "server $(hostname)"
    python server.py
else
    echo "client $(hostname)"
    python worker.py $server
fi
