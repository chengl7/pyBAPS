#!/bin/bash

host=$(hostname)
echo "in test $(hostname)"

if [ "$(hostname)" == "$1" ]; then
    echo "server $(hostname)"
else
    echo "client $(hostname)"
fi
